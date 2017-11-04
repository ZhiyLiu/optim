#ifndef MOVABLE_FL_WINDOW_H
#define MOVABLE_FL_WINDOW_H

/*  Class movable_Fl_Window.

    This class exists because of problems in fluid and FLTK.
    Fluid will write .fl files with calls to the five parameter
    constructor of class Fl_Window, which has arguments of x, y,
    w, h, and title.  However, when fluid generates code, it
    always drops the x and y arguments.  FLTK does not permit
    windows constructed by the three parameter constructor to
    be repositioned.  Windows defined by use of this class
    can, however, be moved, because this class guarantees that
    the five parameter constructor will be used.

    Note that unless a window is repostioned, it will always be
    placed at location (0, 0), if fluid was used to generate the
    constructor call.

    In addition to the basic capability of providing for window
    relocation, this class maintains a list of the windows
    constructed in Pablo.  As instances of this class are
    constructed or deleted, this list is updated, and at any
    time, function madeWindows() can be called to access this
    list.

    A capability is also provided for determining when FL_SHOW
    or FL_HIDE events occur for any movable_Fl_Window by the
    specification of callback functions.

    In addition, a callback is provided for specifying a handle
    function to process other kinds of events.  Note that this
    should not process FL_SHOW and FL_HIDE events, which will
    already have been processed before the provided function is
    called.  Events that are not processed by one of these
    methods, will be passed on the the Fl_Window class.

    Finally, there is also a capability for making adjustments
    after a resize.  Any function pointer passed to function
    setResizeCallback() will be called after the window is
    resized.
*/

#define MAX_NUM_WINDOWS    30


class Fl_Window;


class movable_Fl_Window : public Fl_Window
{

  public:

    movable_Fl_Window(int x, int y, int w, int h, const char *title = 0);
    movable_Fl_Window(int w, int h, const char *title = 0);
    virtual ~movable_Fl_Window();

    void setShowCallback(void (* cb)());
    void setHideCallback(void (* cb)());

    void setResizeCallback(void (* cb)(int, int, int, int));

    void setHandleCallback(int (* cb)(int));

    int handle(int event);
    void resize(int x, int y, int w, int h);

    Fl_Window ** madeWindows(int & len) {
        len = listLen;
        return windowList;
    }

  private:

    static void (* showCallback)();
    static void (* hideCallback)();
    static void (* resizeCallback)(int, int, int, int);
    static int (* handleCallback)(int);
    static Fl_Window * windowList[MAX_NUM_WINDOWS];
    static int listLen;

};



#endif    /* MOVABLE_FL_WINDOW_H */


