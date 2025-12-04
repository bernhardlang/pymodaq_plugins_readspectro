from qtpy.QtWidgets import QMainWindow
from qtpy.QtCore import Qt

class MainWindow(QMainWindow):
    """Extends QMainWindow to catch close events.

    Before closing the application, it may be necessary to perform shutdown
    operations on some hardware. This can be done by registering an appropriate
    callback function which will be executed when a close event is sent to the
    main window. When that callback function returns something else than None,
    the close event is ignored. Othervise it is forwarded to QMainWindow's
    closeEvent method which eventually closes the window.

    Closing the window can be inhibitted to prevent an inappropriate or
    accidental application shutdown. For instance, tearing down hardware
    controllers upon reception of a close event while some view is still
    grabbing data may lead to invalid hardware operations.
    """

    def __init__(self):
        QMainWindow.__init__(self)
        self.close_enabled = True
        self.shutdown_callback = None

    def set_shutdown_callback(self, callback):
        if callback is not None and not callable(callback):
            raise RuntimeError("shutdown callback has type '%s' which is not callable" % type(callback))
        self.shutdown_callback = callback

    # method needed because Python lambdas cannot assign values to attributes
    def enable_close(self, enable):
        self.close_enabled = enable

    def disable_close(self, disable):
        self.enable_close(not disable)

    def closeEvent(self, event):
        if self.close_enabled:
            if self.shutdown_callback is not None:
                close = self.shutdown_callback() is None
            else:
                close = True
            if close:
                QMainWindow.closeEvent(self, event)
                return

        event.ignore() # to be checked if that is needed or even wanted
