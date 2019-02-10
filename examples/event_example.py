#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: event_example
# Author: Glen Langston
# Description: c++ block test
# Generated: Fri Feb  8 13:31:17 2019
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from PyQt4.QtCore import QObject, pyqtSlot
from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import qtgui
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import numpy as np
import radio_astro
import sip
import sys


class event_example(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "event_example")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("event_example")
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "event_example")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())

        ##################################################
        # Variables
        ##################################################
        self.vec_length = vec_length = 512
        self.samp_rate = samp_rate = 10e6
        self.nt = nt = 512
        self.ndm = ndm = 100
        self.freq = freq = 705e6
        self.fftsize = fftsize = 1024
        self.Mode = Mode = 2
        self.Bandwidth = Bandwidth = 2e6

        ##################################################
        # Blocks
        ##################################################
        self._Mode_options = (1, 2, )
        self._Mode_labels = ('Monitor', 'Detect', )
        self._Mode_tool_bar = Qt.QToolBar(self)
        self._Mode_tool_bar.addWidget(Qt.QLabel('Mode'+": "))
        self._Mode_combo_box = Qt.QComboBox()
        self._Mode_tool_bar.addWidget(self._Mode_combo_box)
        for label in self._Mode_labels: self._Mode_combo_box.addItem(label)
        self._Mode_callback = lambda i: Qt.QMetaObject.invokeMethod(self._Mode_combo_box, "setCurrentIndex", Qt.Q_ARG("int", self._Mode_options.index(i)))
        self._Mode_callback(self.Mode)
        self._Mode_combo_box.currentIndexChanged.connect(
        	lambda i: self.set_Mode(self._Mode_options[i]))
        self.top_layout.addWidget(self._Mode_tool_bar)
        self.radio_astro_event_0 = radio_astro.event(fftsize, 3.6, Bandwidth, .001, Mode)
        self.qtgui_time_sink_x_0 = qtgui.time_sink_c(
        	fftsize, #size
        	Bandwidth, #samp_rate
        	"", #name
        	1 #number of inputs
        )
        self.qtgui_time_sink_x_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0.set_y_axis(-1, 1)
        
        self.qtgui_time_sink_x_0.set_y_label('Amplitude', "")
        
        self.qtgui_time_sink_x_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0.enable_grid(False)
        self.qtgui_time_sink_x_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0.enable_control_panel(False)
        
        if not True:
          self.qtgui_time_sink_x_0.disable_legend()
        
        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        
        for i in xrange(2*1):
            if len(labels[i]) == 0:
                if(i % 2 == 0):
                    self.qtgui_time_sink_x_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0.set_line_alpha(i, alphas[i])
        
        self._qtgui_time_sink_x_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_win)
        self.blocks_vector_to_stream_0 = blocks.vector_to_stream(gr.sizeof_gr_complex*1, fftsize)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*fftsize, Bandwidth/fftsize,True)
        self.blocks_stream_to_vector_0_0 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, fftsize)
        self.analog_fastnoise_source_x_0 = analog.fastnoise_source_c(analog.GR_GAUSSIAN, 1, 0, 8192)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_fastnoise_source_x_0, 0), (self.blocks_stream_to_vector_0_0, 0))    
        self.connect((self.blocks_stream_to_vector_0_0, 0), (self.radio_astro_event_0, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.blocks_vector_to_stream_0, 0))    
        self.connect((self.blocks_vector_to_stream_0, 0), (self.qtgui_time_sink_x_0, 0))    
        self.connect((self.radio_astro_event_0, 0), (self.blocks_throttle_0, 0))    

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "event_example")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_vec_length(self):
        return self.vec_length

    def set_vec_length(self, vec_length):
        self.vec_length = vec_length

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate

    def get_nt(self):
        return self.nt

    def set_nt(self, nt):
        self.nt = nt

    def get_ndm(self):
        return self.ndm

    def set_ndm(self, ndm):
        self.ndm = ndm

    def get_freq(self):
        return self.freq

    def set_freq(self, freq):
        self.freq = freq

    def get_fftsize(self):
        return self.fftsize

    def set_fftsize(self, fftsize):
        self.fftsize = fftsize
        self.blocks_throttle_0.set_sample_rate(self.Bandwidth/self.fftsize)

    def get_Mode(self):
        return self.Mode

    def set_Mode(self, Mode):
        self.Mode = Mode
        self._Mode_callback(self.Mode)

    def get_Bandwidth(self):
        return self.Bandwidth

    def set_Bandwidth(self, Bandwidth):
        self.Bandwidth = Bandwidth
        self.qtgui_time_sink_x_0.set_samp_rate(self.Bandwidth)
        self.blocks_throttle_0.set_sample_rate(self.Bandwidth/self.fftsize)


def main(top_block_cls=event_example, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
