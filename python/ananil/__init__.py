
import sys
from soma.qt_gui import ipkernel_tools
import anatomist.direct.api as ana
from .ananil import ContrastPanel

def main():
    # test dataset: /neurospin/unicog/protocols/IRMf/Dighiero_Dehaene_7TLinguistics_2023/

    ipkernel_tools.before_start_ipkernel()

    a = ana.Anatomist(*sys.argv)
    # exits abruptly when closing the control window. Works around the callback
    # problem: when the event loop runs through ipython kernel,
    # QApplication.aboutToQuit signals are never executed, thus cannot actually
    # quit.
    a.setExitOnQuit(True)

    cp = ContrastPanel()
    cp.show()

    ipkernel_tools.start_ipkernel_qt_engine()
