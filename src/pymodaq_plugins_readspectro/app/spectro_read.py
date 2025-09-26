import numpy as np
import csv, time
from qtpy.QtCore import QByteArray, QSettings, QTimer
from qtpy.QtGui import QKeySequence
from qtpy.QtWidgets import QMainWindow, QWidget, QApplication, QProgressBar, \
    QFileDialog
from pyqtgraph import GraphicsLayoutWidget, PlotDataItem, FillBetweenItem
from pyqtgraph import PlotItem, PlotDataItem, ViewBox
from pyqtgraph import GraphicsWidget, PlotWidget
from pymodaq.control_modules.daq_viewer import DAQ_Viewer
from pymodaq.utils.data import DataToExport, DataFromPlugins
from pymodaq_gui.utils.custom_app import CustomApp
from pymodaq_gui.plotting.data_viewers.viewer1D import Viewer1D
from pymodaq_gui.utils.dock import DockArea, Dock
from pymodaq_gui.utils.main_window import MainWindow


RAW             = 0
WITH_BACKGROUND = 1
ABSORPTION      = 2

CONTINUOUS = 0
LINEAR     = 1
LIN_LOG    = 2


class SpectroApp(CustomApp):

    measurement_modes = { 'Raw': RAW, 'Background Subtracted': WITH_BACKGROUND,
                          'Absorption': ABSORPTION }
    scan_modes = { 'Continuous': CONTINUOUS, 'Linear Sequential': LINEAR,
                   'Lin-log Sequential': LIN_LOG }
    scan_params = ['linear_step', 'scan_end', 'log_start', 'points_per_decade' ]
    visible_params = { CONTINUOUS: [], LINEAR: ['linear_step', 'scan_end'],
                       LIN_LOG: scan_params }

    params = [{'name': 'integration_time', 'title': 'Integration Time [sec]',
               'type': 'float', 'min': 0.001, 'max': 100, 'value': 0.05,
               'tip': 'Integration time in seconds'},
              {'name': 'averaging', 'title': 'Averaging',
               'type': 'int', 'min': 1, 'max': 1000, 'value': 1,
               'tip': 'Software Averaging'},
              {'name': 'measurement_mode', 'title': 'Measurement Mode',
               'type': 'list', 'limits': list(measurement_modes.keys()),
               'tip': 'Measurement Mode'},
              {'name': 'scan_mode', 'title': 'Scan Mode', 'type': 'list',
               'limits': list(scan_modes.keys()), 'tip': 'Scan Mode'},
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

    def __init__(self, parent: DockArea, plugin):
        super().__init__(parent)

        self.plugin = plugin
        self.setup_ui()

        # keep screen geometry between runs, should be integrated into
        # PyMoDAQ settings. Is anyway kind of messy because Qt and
        # pyqtgraph don't handle the matter very consistently.
        settings = QSettings("chiphy", "spectro-read")
        geometry = settings.value("geometry", QByteArray())
        self.mainwindow.restoreGeometry(geometry)
        state = settings.value("dockarea", None)
        if state is not None:
            try:
                self.dockarea.restoreState(state)
            except: # pyqtgraph's state restoring is not very fail safe
                # erease inconsistent settings in case pyqtgraph trips
                settings.setValue("dockarea", None)

        # Retrieve spacing of first column in case the user has made it fully
        # visible in a previous run of the program.
        header = settings.value("settings-header-0", None)
        if header is not None:
            self._settings_tree.widget.header().resizeSection(0, int(header))

        self.measurement_mode = RAW
        self.have_background = False
        self.have_reference = False
        self.scan_mode = CONTINUOUS
        self.acquiring = False
        self.adjust_actions()
        self.adjust_parameters()

    def setup_docks(self):
        # left column: essential parameters at top, small plots for dark and
        # reference signals
        # main area for current data

        # top left, essential parameters
        self.docks['settings'] = Dock('Application Settings')
        self.dockarea.addDock(self.docks['settings'])
        self.docks['settings'].addWidget(self.settings_tree)

        # main area with spectrum plot
        spectrum_dock = Dock('Data')
        self.docks['spectrum'] = \
            self.dockarea.addDock(spectrum_dock, "right", self.docks['settings'])
        spectrum_widget = QWidget()
        self.spectrum_viewer = Viewer1D(spectrum_widget)
        spectrum_dock.addWidget(spectrum_widget)

        # progress bar for time sequence
        progress_dock = Dock("Progress")
        self.docks['progress'] = \
            self.dockarea.addDock(progress_dock, 'bottom',
                                  self.docks['settings'])
        progress_widget = QProgressBar()
        progress_dock.addWidget(progress_widget)

        # plot for raw spectrum data and reference 
        raw_data_dock = Dock('Raw Data')
        self.docks['raw-data'] = \
            self.dockarea.addDock(raw_data_dock, "bottom",
                                  self.docks['progress'])
        raw_data_widget = QWidget()
        self.raw_data_viewer = Viewer1D(raw_data_widget)
        self.raw_data_viewer.toolbar.hide()

        raw_data_dock.addWidget(raw_data_widget)

        # plot for background 
        background_dock = Dock('Background')
        self.docks['background'] = \
            self.dockarea.addDock(background_dock, "bottom",
                                  self.docks['raw-data'])
        background_widget = QWidget()
        self.background_viewer = Viewer1D(background_widget)
        background_dock.addWidget(background_widget)
        self.background_viewer.toolbar.hide()

        # separate window with raw detector data
        self.daq_viewer_area = DockArea()
        self.detector = \
            DAQ_Viewer(self.daq_viewer_area, title=self.plugin, init_h5=False)
        self.detector.daq_type = 'DAQ1D'
        self.detector.detector = self.plugin
        self.detector.init_hardware()

    def setup_actions(self):
        self.add_action('acquire', 'Acquire', 'spectrumAnalyzer',
                        "Acquire", checkable=False, toolbar=self.toolbar)
        self.add_action('background', 'Take Background', 'camera',
                        "Take Background", checkable=False, toolbar=self.toolbar)
        self.add_action('reference', 'Take Reference', 'camera',
                        "Take Reference", checkable=False, toolbar=self.toolbar)
        self.add_action('save', 'Save', 'SaveAs', "Save current data",
                        checkable=False, toolbar=self.toolbar)        
        self.add_action('show', 'Show/hide', 'read2', "Show Hide DAQViewer",
                        checkable=True, toolbar=self.toolbar)

    def connect_things(self):
        self.quit_action.triggered.connect(self.mainwindow.close)
        self.connect_action('save', self.save_current_data)
        self.connect_action('show', self.show_detector)
        self.connect_action('acquire', self.start_acquiring)
        self.connect_action('background', self.take_background)
        self.connect_action('reference', self.take_reference)
        self.detector.grab_done_signal.connect(self.show_data)

    def setup_menu(self):
        file_menu = self.mainwindow.menuBar().addMenu('File')
        self.affect_to('save', file_menu)

        file_menu.addSeparator()
        self.quit_action = file_menu.addAction("Quit", QKeySequence('Ctrl+Q'))

    def value_changed(self, param):
        # <<-- param widget should be readonly during measurement
        if param.name() == "integration_time":
            self.detector.settings.child('detector_settings',
                                         'integration_time') \
                                  .setValue(param.value() * 1000)
            # background and reference should be measurement with the same i.t.
            if self.measurement_mode in [WITH_BACKGROUND, ABSORPTION]:
                self.detector.stop()
            self.have_background = False
            self.have_reference = False
            self.adjust_actions()

        if param.name() == "averaging":
            self.detector.settings.child('main_settings', 'Naverage') \
                                         .setValue(param.value())
            self.have_background = False
            self.have_reference = False
            self.adjust_actions()

        elif param.name() == "measurement_mode":
            self.measurement_mode = self.measurement_modes[param.value()]

        elif param.name() == "scan_mode":
            self.detector.stop()
            self.scan_mode = self.scan_modes[param.value()]
            self.adjust_parameters()

        self.adjust_operation()
        self.adjust_actions()

    def adjust_operation(self):
        """Stop acquisition if background / reference is missing but needed"""
        if self.measurement_mode >= WITH_BACKGROUND:
            if not self.have_background:
                self.detector.stop()
            elif self.measurement_mode == ABSORPTION and not self.have_reference:
                self.detector.stop()

    def adjust_actions(self):
        """Disable actions which need other actions to be performed first.
        A reference can only be taken when a background has been measured.
        Acquisition in absorption mode needs a reference (and therefore also
        a background).
        """
        if self.measurement_mode == RAW:
            self._actions["acquire"].setEnabled(True)
            self._actions["background"].setEnabled(False)
            self._actions["reference"].setEnabled(False)
        if self.measurement_mode == WITH_BACKGROUND:
            self._actions["acquire"].setEnabled(self.have_background)
            self._actions["background"].setEnabled(True)
            self._actions["reference"].setEnabled(False)
        if self.measurement_mode == ABSORPTION:
            self._actions["acquire"].setEnabled(self.have_reference)
            self._actions["background"].setEnabled(True)
            self._actions["reference"].setEnabled(self.have_background)

    def adjust_parameters(self):
        """Hide parameters which are not needed in current measurment mode."""
        for child in self.scan_params:
            if child in self.visible_params[self.scan_mode]:
                self.settings.child(child).show()
            else:
                self.settings.child(child).hide()

    def show_detector(self, status):
        self.daq_viewer_area.setVisible(status)

    def show_data(self, data: DataToExport):
        """Display incoming data.

        In linear or logarithmic scan mode: schedule next acquisition
        afterwards.
        ----------
        data: DataToExport
        """
        data1D = data.get_data_from_dim('Data1D')
        signal = data1D[0]
        try:
            ava_time_stamp = data.get_data_from_name('timestamp')[0][0] / 100
        except:
            ava_time_stamp = 0
        system_time_stamp = time.time_ns() * 1e-6

        if self.measurement_mode != RAW:
            if not hasattr(self, 'background'):
                return
            signal[0] -= self.background
            signal.labels = ['signal-background']

        if self.measurement_mode == ABSORPTION:
            if not hasattr(self, 'reference'):
                return
            valid_mask = np.logical_and(signal[0] > 0, self.reference_valid_mask)
            absorption = np.where(valid_mask,
                                  -np.log10(signal[0] / self.reference), 0)
            dfp = DataFromPlugins(name='absorption', data=[absorption],
                              dim='Data1D', labels=['absorption'])
            self.spectrum_viewer.show_data(dfp)
            dfp = DataFromPlugins(name='raw',
                                  data=[signal[0], self.reference],
                                  dim='Data1D', labels=['signal', 'reference'])
            self.raw_data_viewer.show_data(dfp)
            self.current_data = absorption
        else:
            self.current_data = signal[0]
            self.spectrum_viewer.show_data(signal)

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

        if self.acquiring: # rather stop it
            self.acquiring = False
            self.detector.stop_grab()
            return

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

    def take_background(self):
        """Grab one background spectrum."""

        if self.detector.controller.has_dark_shutter:
            self.detector.controller.open_dark_shutter(False)

        n_average = \
            self.detector.settings.child('main_settings', 'Naverage').value()
        self.background,timestamp = self.detector.controller.grab_spectrum()
        for _ in range(1, n_average):
            data,timestamp = self.detector.controller.grab_spectrum()
            self.background += data
        self.background /= n_average
        self.have_background = True
        self.adjust_actions()
        dfp = DataFromPlugins(name='Avantes', data=self.background, dim='Data1D',
                              labels=['background'])
        self.spectrum_viewer.show_data(dfp)
        self.background_viewer.show_data(dfp)
        if self.detector.controller.has_dark_shutter:
            self.detector.controller.open_dark_shutter(True)

    def take_reference(self):
        """Grab one reference spectrum."""

        if self.detector.controller.has_reference_switch:
            self.detector.controller.set_reference_switch(True)

        n_average = \
            self.detector.settings.child('main_settings', 'Naverage').value()

        self.reference = np.zeros(len(self.background))
        self.reference_valid_mask = \
            np.full(len(self.background), True, dtype=bool)
        for _ in range(n_average):
            data,timestamp = self.detector.controller.grab_spectrum()
            data -= self.background
            self.reference_valid_mask = \
                np.logical_and(data > 0, self.reference_valid_mask)
            self.reference += data

        self.reference /= n_average
        self.have_reference = True
        self.adjust_actions()
        dfp = DataFromPlugins(name='Avantes', data=self.reference, dim='Data1D',
                              labels=['data'])
        self.spectrum_viewer.show_data(dfp)
        dfp = DataFromPlugins(name='Avantes', data=self.reference, dim='Data1D',
                              labels=['reference'])
        self.raw_data_viewer.show_data(dfp)

        if self.detector.controller.has_reference_switch:
            self.detector.controller.set_reference_switch(False)


    def save_current_data(self):
        """Save dat currently displayed on the main plot."""
        if not hasattr(self, 'current_data'):
            return
        result = QFileDialog.getSaveFileName(caption="Save Data", dir=".",
                                             filter="*.csv")
        if result is None or not len(result[0]):
            return
        wavelengths = self.detector.controller.wavelengths
        with open(result[0], "wt") as csv_file:
            writer = csv.writer(csv_file, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for i,wl in enumerate(wavelengths):
                writer.writerow([wl, self.current_data[i]])

    def quit_function(self):
        self.clean_up()
        self.mainwindow.close()

    def clean_up(self):
        self.detector.quit_fun()
        QApplication.processEvents()
        settings = QSettings("chiphy", "spectro-read")
        settings.setValue("geometry", self.mainwindow.saveGeometry())
        settings.setValue("dockarea", self.dockarea.saveState())
        settings.setValue("settings-header-0",
                          self._settings_tree.widget.header().sectionSize(0))


def main():
    import sys
    from pymodaq_gui.utils.utils import mkQApp
    from qtpy.QtCore import pyqtRemoveInputHook

    if len(sys.argv) > 1:
        if sys.argv[1] == '--simulate':
            plugin = "MockSpectro"
        elif len(sys.argv) > 2 and sys.argv[1] == '--plugin':
            plugin=sys.argv[2]
            del sys.argv[2]
            del sys.argv[1]
        else:
            raise RuntimeError("command line argument error")
    else:
        plugin='Avantes'

    app = mkQApp(plugin)
    pyqtRemoveInputHook() # needed for using pdb inside the qt eventloop

    mainwindow = MainWindow()
    dockarea = DockArea()
    mainwindow.setCentralWidget(dockarea)

    prog = SpectroApp(dockarea, plugin=plugin)
    mainwindow.application = prog # not very clean, could be done by event filter
    mainwindow.show()

    app.exec()


if __name__ == '__main__':
    main()
