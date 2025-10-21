import numpy as np
import time
from pymodaq.utils.daq_utils import ThreadCommand
from pymodaq.utils.data import DataFromPlugins, Axis, DataToExport, DataRaw
from pymodaq.control_modules.viewer_utility_classes import DAQ_Viewer_base, \
    comon_parameters, main
from pymodaq.utils.parameter import Parameter
from pymodaq_plugins_avantes.hardware.AvaSpec_ULS2048CL_EVO_Controller \
    import AvantesController, AvantesSimuController


class DAQ_1DViewer_MockSpectro(DAQ_Viewer_base):
    """Simulated Spectrometer Instrument plugin class for a 1D viewer.

    Besides testing spectrometer applications, this plugin permits to simulate
    real detectors with realistic noise settings.
    """

    # no controller, all simulation is done within this class
    controller_type = None

    params = comon_parameters+[
        {'title': 'Integration time [sec]', 'name': 'integration_time',
         'type': 'float', 'min': 0.001, 'value': 1,
         'tip': 'Integration time in seconds' },
        {'title': 'Number of pixels', 'name': 'n_pixels', 'type': 'int',
         'min': 1, 'value': 500 },
        {'title': 'Readout noise [LSB]', 'name': 'readout_noise',
         'type': 'float', 'min': 0, 'value': 4,
         'tip': 'Acquisition noise in LSB' },
        {'title': 'Dark level [LSB/µs]:', 'name': 'dark_level',
         'type': 'float', 'min': 0, 'value': 2e3,
         'tip': 'Dark signal per second in LSB' },
        {'title': 'Light level [LSB]', 'name': 'light_level',
         'type': 'float', 'min': 0, 'value': 32e3,
         'tip': 'Signal per second in LSB' },
        {'title': 'Conversion [PE/LSB]', 'name': 'pe_per_lsb',
         'type': 'float', 'min': 1, 'value': 18.3,
         'tip': 'Photo electrons per LSB' },
        {'title': 'ADC bits', 'name': 'adc_bits',
         'type': 'int', 'min': 10, 'max': 32, 'value': 16,
         'tip': 'Number of ADC bits' },
        {'title': 'Wavelengths from:', 'name': 'wl_from',
         'type': 'float', 'min': 250, 'max': 450, 'value': 300 },
        {'title': 'Wavelengths to:', 'name': 'wl_to',
         'type': 'float', 'min': 450, 'max': 1500, 'value': 900 },
        {'title': 'Absorption:', 'name': 'absorption',
         'type': 'float', 'min': 0, 'value': 0.3,
         'tip': 'Simulated optical density' },
    ]

    def ini_attributes(self):
        self.x_axis = None
        # virtual shutter support
        self.dark_shutter_open = True
        self.reference_switch = True # bypass switch around sample cell

    def commit_settings(self, param: Parameter):
        pass

    def ini_detector(self, controller=None):
        """Detector communication initialization

        Parameters
        ----------
        controller: (object)
            custom object of a PyMoDAQ plugin (Slave case). None if only 
            one actuator/detector by controller
            (Master case)

        Returns
        -------
        info: str
        initialized: bool
            False if initialization failed otherwise True
        """

        self.ini_detector_init(slave_controller=controller)

        if self.is_master:
            self.controller = self # we'll provide the necessary ourselves
            n_pixels = self.settings.child("n_pixels").value()
            self.wavelengths = \
                np.linspace(self.settings.child("wl_from").value(),
                            self.settings.child("wl_to").value(), n_pixels)
            pixels = np.linspace(0, n_pixels - 1, n_pixels)
            self.spectrum = \
                np.exp(-((pixels - n_pixels / 2) / (n_pixels / 3))**4)
            self.absorption = \
                self.settings.child("absorption").value() \
                * np.exp(-((pixels - n_pixels / 4) / (n_pixels / 8))**2)

            self.x_axis = Axis(label='Wavelength', units='nm',
                                data=self.wavelengths, index=0)
            #dfp = DataFromPlugins(name='MockSpectro',
            #                      data=[np.zeros(len(self.wavelengths))],
            #                      dim='Data1D', axes=[self.x_axis],
            #                      labels=['MockSpectro-Signal'])
            #self.dte_signal_temp.emit(DataToExport(name='MockSpectro',
            #                                       data=[dfp]))

        return "MockSpectro initialised", True

    def close(self):
        """Terminate the communication protocol"""

        self.controller.close_communication()
        initialized = False

    def grab_data(self, Naverage=1, **kwargs):
        """Start grabbing from the detector
        Use a synchrone acquisition (blocking function)

        Parameters
        ----------
        Naverage: int
            Number of hardware averaging.
        """

        spectrum, time_stamp = self.simulate_spectrum()

#        dwa0D_timestamp = \
#            DataRaw('timestamp', units='dimensionless',
#                    data=np.array([timestamp]))

        dfp_spectrum = \
            DataFromPlugins(name='MockSpectro', data=spectrum, dim='Data1D',
                            labels=['spectrum'], axes=[self.x_axis])
        dfp_time_stamp = \
            DataFromPlugins(name='TimeStamp', data=time_stamp, dim='Data0D',
                            labels=['time stamp'])
        
        self.dte_signal.emit(DataToExport(name='spectrum',
                                          data=[dfp_spectrum, dfp_time_stamp]))

    def simulate_spectrum(self):
        n_pixels = self.settings.child("n_pixels").value()
        integration = self.settings.child("integration_time").value()
        time.sleep(integration)
        data = \
            np.random.normal(loc=self.settings.child("dark_level").value()
                             * integration,
                             scale=self.settings.child("readout_noise").value(),
                             size=n_pixels)
        if self.dark_shutter_open:
            pe_per_lsb = self.settings.child("pe_per_lsb").value()
            light = self.spectrum * self.settings.child("light_level").value() \
                * integration * pe_per_lsb
            if not self.reference_switch:
                light *= np.power(10, -self.absorption)
            data += (data + np.random.poisson(light)) / pe_per_lsb
        max_adc = 1 << self.settings.child("adc_bits").value()
        data = np.where(data < max_adc, np.floor(data), max_adc)

        return data, time.time()

    def grab_spectrum(self):
        return self.simulate_spectrum()

    def stop(self):
        pass

    # dummy functions which avoid writing a separate dummy controller which
    # controls ... nothing
    def open_communication(self):
        pass

    def close_communication(self):
        pass

    def has_dark_shutter(self):
        return True

    def open_dark_shutter(self, open_shutter):
        self.dark_shutter_open = open_shutter

    def has_reference_switch(self):
        return True

    def set_reference_switch(self, set_switch):
        self.reference_switch = set_switch


if __name__ == '__main__':
    main(__file__)
