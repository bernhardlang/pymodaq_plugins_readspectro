import numpy as np
import csv, time
from copy import deepcopy
from os import path
from qtpy.QtCore import QByteArray, QSettings, QTimer
from qtpy.QtGui import QKeySequence, QPixmap
from qtpy.QtWidgets import QMainWindow, QWidget, QApplication, QProgressBar, \
    QFileDialog, QSplashScreen
from pyqtgraph import GraphicsLayoutWidget, PlotDataItem, FillBetweenItem
from pyqtgraph import PlotItem, PlotDataItem, ViewBox
from pyqtgraph.dockarea import DockLabel
from pyqtgraph import GraphicsWidget, PlotWidget
from pymodaq.utils.data import DataToExport, DataFromPlugins
from pymodaq_gui.utils.custom_app import CustomApp
from pymodaq_gui.plotting.data_viewers.viewer1D import Viewer1D
from pymodaq_gui.utils.dock import DockArea, Dock
from pymodaq_plugins_readspectro.gui_utils import MainWindow, Custom_DAQ_Viewer


RAW             = 0
WITH_BACKGROUND = 1
ABSORPTION      = 2

class SpectroApp(CustomApp):

    measurement_modes = { 'Raw': RAW, 'Background Subtracted': WITH_BACKGROUND,
                          'Absorption': ABSORPTION }

    params = [{'name': 'integration_time', 'title': 'Integration Time [ms]',
               'type': 'float', 'min': 0.001, 'max': 10000, 'value': 50,
               'tip': 'Integration time in seconds'},
              {'name': 'averaging', 'title': 'Averaging',
               'type': 'int', 'min': 1, 'max': 1000, 'value': 10,
               'tip': 'Software Averaging'},
              {'name': 'pymo_averaging', 'title': 'PyMoDAQ Averaging',
               'type': 'int', 'min': 1, 'max': 1000, 'value': 1,
               'tip': 'Background Software Averaging'},
              {'name': 'back_averaging', 'title': 'Background Averaging',
               'type': 'int', 'min': 1, 'max': 1000, 'value': 100,
               'tip': 'Background Software Averaging'},
              {'name': 'ref_averaging', 'title': 'Reference Averaging',
               'type': 'int', 'min': 1, 'max': 1000, 'value': 100,
               'tip': 'Reference Software Averaging'},
              {'name': 'measurement_mode', 'title': 'Measurement Mode',
               'type': 'list', 'limits': list(measurement_modes.keys()),
               'tip': 'Measurement Mode'},
              ]

    app_name = "absorption"

    def __init__(self, parent: DockArea, plugin, main_window=None):
        super().__init__(parent)

        self.plugin = plugin
        self.setup_ui()

        # keep screen geometry between runs, should be integrated into
        # PyMoDAQ settings. Is anyway kind of messy because Qt and
        # pyqtgraph don't handle the matter very consistently.
        self.qt_settings = QSettings("chiphy", self.app_name)
        geometry = self.qt_settings.value("geometry", QByteArray())
        self.mainwindow.restoreGeometry(geometry)
        state = self.qt_settings.value("dockarea", None)
        if state is not None:
            try:
                self.dockarea.restoreState(state)
            except: # pyqtgraph's state restoring is not very fail safe
                # erease inconsistent settings in case pyqtgraph trips
                self.qt_settings.setValue("dockarea", None)

        # Retrieve spacing of first column in case the user has made it fully
        # visible in a previous run of the program.
        header = self.qt_settings.value("settings-header-0", None)
        if header is not None:
            self._settings_tree.widget.header().resizeSection(0, int(header))

        self.measurement_mode = RAW
        self.have_background = False
        self.have_reference = False
        self.acquiring = False
#        self.n_average = 1
#        self.n_back = 1
#        self.n_ref = 1
        self.adjust_actions()
#        self.delimiter = '\t'
        self.delimiter = ','

        if main_window is not None:
            main_window.set_shutdown_callback(self.quit_function)

    def setup_docks(self):
        # left column: essential parameters at top, small plots for dark and
        # reference signals
        # main area for current data

        # top left, essential parameters
        self.docks['settings'] = Dock('Application Settings')
        self.dockarea.addDock(self.docks['settings'])
        self.docks['settings'].addWidget(self.settings_tree)

        # main area with spectrum plot
        self.spectrum_label = DockLabel("Raw Data")
        spectrum_dock = Dock('Data', label=self.spectrum_label)
        self.docks['spectrum'] = \
            self.dockarea.addDock(spectrum_dock, "right",
                                  self.docks['settings'])
        spectrum_widget = QWidget()
        self.spectrum_viewer = Viewer1D(spectrum_widget)
        self.spectrum_viewer.toolbar.hide()
        spectrum_dock.addWidget(spectrum_widget)

        # plot for raw spectrum data and reference 
        raw_data_dock = Dock('Raw Data')
        self.docks['raw-data'] = \
            self.dockarea.addDock(raw_data_dock, "bottom",
                                  self.docks['settings'])
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
            Custom_DAQ_Viewer(self.daq_viewer_area, title=self.plugin,
                              no_continuous_saving=True)
        self.detector.daq_type = 'DAQ1D'
        self.detector.detector = self.plugin
        self.detector.init_hardware()

        # is this the right place / way?
        self.detector.settings.child('detector_settings', 'integration_time')\
            .setValue(self.settings.child('integration_time').value())

        #self.detector.grab_status.connect(self.mainwindow.disable_close)
        
    def setup_actions(self):
        self.add_action('acquire', 'Acquire', 'run2',
                        "Acquire", checkable=False, toolbar=self.toolbar)
        self.add_action('stop', 'Stop', 'stop',
                        "Stop", checkable=False, toolbar=self.toolbar)
        self.add_action('background', 'Take Background', './half-moon.png',
                        "Take Background", checkable=False, toolbar=self.toolbar)
        self.add_action('reference', 'Take Reference', './sun.png',
                        "Take Reference", checkable=False, toolbar=self.toolbar)
        self.add_action('save', 'Save', 'SaveAs', "Save current data",
                        checkable=False, toolbar=self.toolbar)        
        self.add_action('show', 'Show/hide', 'read2', "Show Hide DAQViewer",
                        checkable=True, toolbar=self.toolbar)
        self._actions["stop"].setEnabled(False)

    def connect_things(self):
        self.quit_action.triggered.connect(self.mainwindow.close)
        self.connect_action('save', self.save_current_data)
        self.connect_action('show', self.show_detector)
        self.connect_action('acquire', self.start_acquiring)
        self.connect_action('stop', self.stop_acquiring)
        self.connect_action('background', self.take_background)
        self.connect_action('reference', self.take_reference)
        self.detector.grab_done_signal.connect(self.take_data)

    def setup_menu(self):
        file_menu = self.mainwindow.menuBar().addMenu('File')
        self.affect_to('save', file_menu)

        file_menu.addSeparator()
        self.quit_action = file_menu.addAction("Quit", QKeySequence('Ctrl+Q'))

    def value_changed(self, param):
        if param.name() == "integration_time":
            self.detector.settings.child('detector_settings',
                                         'integration_time') \
                                  .setValue(param.value())
            # background and reference should be measurement with the same i.t.
            if self.measurement_mode in [WITH_BACKGROUND, ABSORPTION]:
                self.detector.stop()
            self.have_background = False
            self.have_reference = False
            self.adjust_actions()

        if param.name() == "averaging":
            self.average = param.value()
        if param.name() == "pymo_averaging":
            try:
                self.detector.settings.child('main_settings', 'Naverage') \
                                      .setValue(param.value())
            except:
                import pdb
                pdb.set_trace()
        elif param.name() == "back_averaging":
            self.background_average = param.value()
        elif param.name() == "measurement_mode":
            self.measurement_mode = self.measurement_modes[param.value()]

        self.adjust_operation()
        self.adjust_actions()

    def adjust_operation(self):
        """Stop acquisition if background / reference is missing but needed"""
        if self.measurement_mode < WITH_BACKGROUND:
            dock_title = "Raw Data"
        else:
            dock_title = "Absorption" if self.measurement_mode == ABSORPTION \
                else "Background Subtracted Data"
            if not self.have_background:
                self.detector.stop()
            elif self.measurement_mode == ABSORPTION and not self.have_reference:
                self.detector.stop()

        self.spectrum_label.setText(dock_title)

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

    def show_detector(self, status):
        self.daq_viewer_area.setVisible(status)

    def average_data(self, sum_data, squares_data, samples):
        mean = sum_data / samples
        error = np.sqrt((samples * squares_data - sum_data**2)
                        / (samples**2 * (samples - 1)))
        return mean, error

    def accumulate_data(self, data, n_samples):
        if n_samples:
            self.sum_data += data
            self.squares_data += data**2
        else:
            self.sum_data = data
            self.squares_data = data**2
        return n_samples + 1

    def take_data(self, data: DataToExport):
        data1D = data.get_data_from_dim('Data1D')
        self.n_samples = self.accumulate_data(data[0][0], self.n_samples)
        if self.n_samples <= self.n_average:
            return

        self.mean_raw, self.error_raw = \
            self.average_data(self.sum_data, self.squares_data, self.n_samples)
        self.n_samples = 0

        if self.measurement_mode == RAW:
            self.show_data(self.mean_raw, self.error_raw, 'raw')
            return

        self.mean_signal = self.mean_raw - self.background
        self.error_signal = np.sqrt(self.error_raw**2 + self.error_background**2)
        if self.measurement_mode == WITH_BACKGROUND:
            self.show_data(self.mean_signal, self.error_signal,
                           'signal', self.mean_raw)
        else: # self.measurement_mode == ABSORPTION:
            valid_mask = \
                np.logical_and(self.mean_signal > 0, self.reference_valid_mask)
            self.absorption = \
                np.where(valid_mask,
                         -np.log10(self.mean_signal / self.reference), 0)
            self.error_absorption = \
                1 / np.log(10) \
                * np.sqrt((self.error_raw / self.mean_signal)**2
                          + ((self.error_reference + self.error_background)
                             / self.reference)**2
                          + (1 / self.mean_signal - 1 / self.reference)**2
                            * self.error_background)

            self.show_data(self.absorption, self.error_absorption, 'absorption',
                           self.mean_raw, self.reference)

    def show_data(self, mean, error, name, raw=None, reference=None):
        dfp = DataFromPlugins(name=name, data=[mean, error], dim='Data1D',
                              labels=[name, 'error'])
        self.spectrum_viewer.show_data(dfp)
        if raw is not None:
            data = [raw]
            labels = ['raw signal']
            if reference is not None:
                data.append(reference)
                labels.append('reference')
            dfp = DataFromPlugins(name='raw', data=data, dim='Data1D',
                                  labels=labels)
            self.raw_data_viewer.show_data(dfp)

    def start_acquiring(self):
        """Start acquisition"""

        self.docks['settings'].setEnabled(False)
        self._actions["acquire"].setEnabled(False)
        self._actions["stop"].setEnabled(True)
        self.acquiring = True
        self.n_average = self.settings.child('averaging').value()
        self.n_samples = 0
        self.detector.grab()

    def stop_acquiring(self):
        self.docks['settings'].setEnabled(True)
        self._actions["acquire"].setEnabled(True)
        self._actions["stop"].setEnabled(False)
        self.acquiring = False
        self.detector.stop_grab()

    def take_background(self):
        """Grab one background spectrum."""

        if hasattr(self.detector.controller, "open_dark_shutter"):
            self.detector.controller.open_dark_shutter(False)

        self.n_back = self.settings.child('back_averaging').value()
        for i in range(0, self.n_back):
            background,timestamp = self.detector.controller.grab_spectrum()
            self.accumulate_data(background, i)
        self.background, self.error_background = \
            self.average_data(self.sum_data, self.squares_data, self.n_back)
        self.have_background = True
        self.adjust_actions()
        dfp = DataFromPlugins(name='Avantes',
                              data=[self.background, self.error_background],
                              dim='Data1D', labels=['background', 'error'])
        self.spectrum_viewer.show_data(dfp)
        self.background_viewer.show_data(dfp)
        if hasattr(self.detector.controller, "open_dark_shutter"):
            self.detector.controller.open_dark_shutter(True)

    def take_reference(self):
        """Grab one reference spectrum."""

        if hasattr(self.detector.controller, "set_reference_switch"):
            self.detector.controller.set_reference_switch(True)

        self.n_ref = self.settings.child('ref_averaging').value()

        for i in range(0, self.n_back):
            reference,timestamp = self.detector.controller.grab_spectrum()
            self.accumulate_data(reference, i)
        self.reference, self.error_reference = \
            self.average_data(self.sum_data, self.squares_data, self.n_ref)
        self.reference -= self.background

        self.reference_valid_mask = self.reference > 0
        self.have_reference = True
        self.adjust_actions()
        dfp = DataFromPlugins(name='Avantes',
                              data=[self.reference, self.error_reference],
                              dim='Data1D', labels=['reference', 'error'])
        self.spectrum_viewer.show_data(dfp)
        dfp = DataFromPlugins(name='Avantes',
                              data=[self.reference, self.error_reference],
                              dim='Data1D', labels=['reference', 'error'])
        self.raw_data_viewer.show_data(dfp)

        if hasattr(self.detector.controller, "set_reference_switch"):
            self.detector.controller.set_reference_switch(False)

    def save_current_data(self):
        """Save dat currently displayed on the main plot."""
        directory = self.qt_settings.value('data-dir', None)
        if directory is None:
            directory = "."
        result = QFileDialog.getSaveFileName(caption="Save Data", dir=directory,
                                             filter="*.csv")
        if result is None or not len(result[0]):
            return

        self.qt_settings.setValue('data-dir', path.dirname(result[0]))

        wavelengths = self.detector.controller.wavelengths
        with open(result[0], "wt") as csv_file:
            writer = csv.writer(csv_file, delimiter=self.delimiter,
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
            if self.measurement_mode == RAW or not self.have_background:
                writer.writerow(['wavelength', 'raw data', 'error'])
                for i,wl in enumerate(wavelengths):
                    writer.writerow(['%.1f' % wl, '%.3f' % self.mean_raw[i],
                                    '%.3f' % self.error_raw[i]])
                return

            if self.measurement_mode == WITH_BACKGROUND \
               or not self.have_reference:
                writer.writerow(['wavelength', 'raw data', 'error raw',
                                 'background', 'error background',
                                 'background subtracted', 'error'])
                for i,wl in enumerate(wavelengths):
                    writer.writerow(['%.1f' % wl, '%.3f' % self.mean_raw[i],
                                     '%.1f' % self.error_raw[i],
                                     '%.1f' % self.background[i],
                                     '%.1f' % self.error_background[i],
                                     '%.1f' % self.mean_signal[i],
                                     '%.1f' % self.error_signal[i]])
                return

            # self.measurement_mode == ABSORPTION
            writer.writerow(['wavelength', 'raw data', 'error raw', 'background',
                             'error background', 'reference', 'error reference',
                             'absorption', 'error'])
            for i,wl in enumerate(wavelengths):
                writer.writerow(['%.1f' % wl, '%.3f' % self.mean_raw[i],
                                 '%.3f' % self.error_raw[i],
                                 '%.3f' % self.background[i],
                                 '%.3f' % self.error_background[i],
                                 '%.3f' % self.reference[i],
                                 '%.3f' % self.error_reference[i],
                                 '%.6f' % self.absorption[i],
                                 '%.6f' % self.error_absorption[i]])

    def quit_function(self):
        self.clean_up()
        self.mainwindow.close()

    def clean_up(self):
        self.detector.quit_fun()
        QApplication.processEvents()
        self.qt_settings.setValue("geometry", self.mainwindow.saveGeometry())
        self.qt_settings.setValue("dockarea", self.dockarea.saveState())
        self.qt_settings.setValue("settings-header-0",
                          self._settings_tree.widget.header().sectionSize(0))


def main(prog):
    import sys, time
    from pymodaq_gui.utils.utils import mkQApp
    from qtpy.QtCore import pyqtRemoveInputHook, QTimer
    from qtpy.QtWidgets import QSplashScreen


    app = QApplication(sys.argv)
#    app.processEvents()

    show_splash = True
    plugin = 'Avantes'

    arg_pos = 1
    while arg_pos < len(sys.argv):
        if sys.argv[arg_pos] == '--simulate':
            plugin = "MockSpectro"
        elif sys.argv[arg_pos] == '--plugin' and arg_pos < len(sys.argv) - 1:
            arg_pos += 1
            plugin=sys.argv[arg_pos]
        elif sys.argv[arg_pos] == '--no-splash':
            show_splash = False
        else:
            raise RuntimeError("command line argument error")
        arg_pos += 1

    if show_splash:
        splash_pixmap = QPixmap("splash.png")
        splash = QSplashScreen(splash_pixmap)
        splash.show()

    app = mkQApp(plugin)
    pyqtRemoveInputHook() # needed for using pdb inside the qt eventloop

    mainwindow = MainWindow()
    dockarea = DockArea()
    mainwindow.setCentralWidget(dockarea)

    prog = prog(dockarea, plugin=plugin, main_window=mainwindow)

    if show_splash:
        QTimer.singleShot(2000, splash.close)
        QTimer.singleShot(2000, mainwindow.show)
    else:
        mainwindow.show()

    sys.exit(app.exec())

if __name__ == '__main__':
    main(SpectroApp)
