from qtpy.QtWidgets import QProgressBar, QMainWindow
from pymodaq_gui.utils.dock import DockArea, Dock
from pymodaq.utils.data import DataToExport
from pymodaq_plugins_readspectro.app.absorption import SpectroApp, main


CONTINUOUS = 0
LINEAR     = 1
LIN_LOG    = 2


class TimedSpectroApp(SpectroApp):

    scan_modes = { 'Continuous': CONTINUOUS, 'Linear Sequential': LINEAR,
                   'Lin-log Sequential': LIN_LOG }
    scan_params = ['reference_refresh', 'linear_step', 'scan_end', 'log_start',
                   'points_per_decade' ]
    visible_params = { CONTINUOUS: [],
                       LINEAR: ['reference_refresh', 'linear_step', 'scan_end'],
                       LIN_LOG: scan_params }

    params = SpectroApp.params + \
        [{'name': 'scan_mode', 'title': 'Scan Mode', 'type': 'list',
          'limits': list(scan_modes.keys()), 'tip': 'Scan Mode'},
         {'name': 'reference_refresh', 'title': 'Refresh Reference [sec]',
          'type': 'float', 'min': 1, 'max': 1000, 'value': 10,
          'tip': 'Minimum interval after which the reference is refreshed'},
         {'name': 'linear_step', 'title': 'Linear Scan Step [sec]',
          'type': 'float', 'min': 0.001, 'max': 1000, 'value': 1,
          'tip': 'Scan step in seconds'},
         {'name': 'scan_end', 'title': 'Scan End [sec]', 'value': 100,
          'type': 'float', 'min': 1, 'max': 10000,
          'tip': 'End of measurement in seconds'},
         {'name': 'log_start', 'title': 'Log Start [sec]',
          'type': 'float', 'min': 1, 'max': 10000, 'value': 10,
          'tip': 'Start of logarithmic scale in seconds'},
         {'name': 'points_per_decade', 'title': 'Points per Decade',
          'type': 'float', 'min': 1, 'max': 10000, 'value': 10,
          'tip': 'Logarithmic points per time decade measurement'},
         ]

    app_name = "timed-absorption"
    
    def __init__(self, parent: DockArea, plugin, main_window=None):
        super().__init__(parent, plugin, main_window=main_window)
        self.scan_mode = CONTINUOUS
        self.adjust_parameters()

    def setup_docks(self):
        super().setup_docks()
        # progress bar for time sequence
        progress_dock = Dock("Progress")
        self.docks['progress'] = \
            self.dockarea.addDock(progress_dock, 'bottom',
                                  self.docks['settings'])
        progress_widget = QProgressBar()
        progress_dock.addWidget(progress_widget)
        
    def value_changed(self, param):
        if param.name() == "scan_mode":
            self.detector.stop()
            self.scan_mode = self.scan_modes[param.value()]
            self.adjust_parameters()
        super().value_changed(param)

    def adjust_parameters(self):
        """Hide parameters which are not needed in current measurment mode."""
        for child in self.scan_params:
            if child in self.visible_params[self.scan_mode]:
                self.settings.child(child).show()
            else:
                self.settings.child(child).hide()

    def show_data(self, data: DataToExport):
        try:
            ava_time_stamp = data.get_data_from_name('timestamp')[0][0] / 100
        except:
            ava_time_stamp = 0
        system_time_stamp = time.time_ns() * 1e-6

        super().show_data(data)
        if self.scan_mode == CONTINUOUS or not self.acquiring:
            return

        # scanning mode: store data and reschedule
        if self.system_start_time is None:
            self.ava_start_time = ava_time_stamp
            self.system_start_time = system_time_stamp
        ava_time_stamp -= self.ava_start_time
        system_time_stamp -= self.system_start_time

        if self.measurement_mode == ABSORPTION:
            # convert absorption data to mOD
            self.write_spectrum(int(system_time_stamp + 0.5),
                                int(ava_time_stamp + 0.5),
                                self.current_data * 1000)
        else:
            self.write_spectrum(int(system_time_stamp + 0.5),
                                int(ava_time_stamp + 0.5), self.current_data)

        # advance time index until next measurement is in future
        old_idx = self.current_time_index
        while self.scheduled_measurement_times[self.current_time_index] \
              <= system_time_stamp:
            self.current_time_index += 1
            if self.current_time_index == len(self.scheduled_measurement_times):
                self.data_file.close() # done
                return

        if self.current_time_index == old_idx + 1:
            # schedule next measurement
            time_delay = \
                self.scheduled_measurement_times[self.current_time_index] \
                - system_time_stamp
            QTimer.singleShot(int(time_delay), self.detector.snap)
        else: # missed scheduled point(s): restart measurement immediately
            self.detector.snap()

    def start_acquiring(self):
        """Start acquisition"""

        self._actions["acquire"].setEnabled(False)
        self._actions["stop"].setEnabled(True)
        self.acquiring = True
        if self.scan_mode == CONTINUOUS:
            self.detector.grab() # just go
            return

        result = QFileDialog.getSaveFileName(caption="Save Data", dir=".",
                                             filter="*.csv")
        if result is None:
            return

        # <<-- unclear: where does init of param values happen?
        
        #self.detector.settings.child('detector_settings', 'integration_time') \
        #    .setValue(self._settings_tree.value("integration_time") * 1000)
        
        # determine number of significant digits according to
        # error = sqrt(Naverage) assuming the best case with
        # error(Naverage=1) is 1
        n_average = \
            self.detector.settings.child('main_settings', 'Naverage').value()
        if n_average > 1:
            self.format_string = \
                '\t{val:.%df}' % (int(np.log10(np.sqrt(n_average))) + 1)
        else:
            self.format_string = '\t{val:.0f}'

        # write wavelengths into file
        self.data_file = open(result[0], "wt")
        self.data_file.write('0\t0')
        for wl in self.detector.controller.wavelengths:
            self.data_file.write('\t%.1f' % wl)
        self.data_file.write('\n')

        # background (time delay -1) and reference (time delay -2) if needed
        if self.measurement_mode >= WITH_BACKGROUND:
            self.write_spectrum(-1, 0, self.background)
            if self.measurement_mode == ABSORPTION:
                self.write_spectrum(-2, 0, self.reference)

        # calculate schedule, only linear for the moment, ignore lin-log
        # QTimer counts in milliseconds
        scan_end = self.settings['scan_end'] * 1000
        linear_step = self.settings['linear_step'] * 1000
        n_measurements = int(scan_end / linear_step) + 1
        self.scheduled_measurement_times = \
            np.linspace(0, scan_end, n_measurements) \
            - self.settings['integration_time'] * 1000 - 20

        # and go
        self.current_time_index = 0
        self.system_start_time = None
        self.detector.snap()

    def write_spectrum(self, t1, t2, spectrum):
        """Writes a single spectrum to file.
        The first two columns contain the system time at data retrieval
        and the time stamp returned by the avaspec library, respectively.
        """
        self.data_file.write('%d\t%d'
                             % (int(np.floor(t1 + 0.5)),
                                int(np.floor(t2 + 0.5))))
        if self.measurement_mode == ABSORPTION and t1 != -1:
            # To save space in the data file, absorption data are stored as
            # integer numbers. 1 LSB is 1µOD.
            for value in spectrum:
                self.data_file.write('\t%d' % (value * 1000 + 0.5))
        else:
            for value in spectrum:
                self.data_file.write(self.format_string.format(val = value))
        self.data_file.write('\n')


if __name__ == '__main__':
    main(TimedSpectroApp)
