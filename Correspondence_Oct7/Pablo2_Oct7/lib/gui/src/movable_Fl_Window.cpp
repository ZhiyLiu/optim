
#include <stdlib.h>
#include <FL/Fl_Window.H>
#include "movable_Fl_Window.h"

//#define DEBUG

#ifdef DEBUG
#include <iostream>
#include <assert.h>
using namespace std;
#endif

void (* movable_Fl_Window::hideCallback)() = NULL;
void (* movable_Fl_Window::showCallback)() = NULL;
void (* movable_Fl_Window::resizeCallback)(int, int, int, int) = NULL;
int (* movable_Fl_Window::handleCallback)(int) = NULL;
Fl_Window * movable_Fl_Window::windowList[MAX_NUM_WINDOWS];
int movable_Fl_Window::listLen = 0;


movable_Fl_Window::movable_Fl_Window(int x, int y, int w, int h, const char *title)
    : Fl_Window(x, y, w, h, title)
{
#ifdef DEBUG
    assert(listLen < MAX_NUM_DOUBLE_WINDOWS - 1);
#endif
    windowList[listLen++] = this;
}

movable_Fl_Window::movable_Fl_Window(int w, int h, const char *title)
    : Fl_Window(0, 0, w, h, title)
{
#ifdef DEBUG
    assert(listLen < MAX_NUM_DOUBLE_WINDOWS - 1);
#endif
    windowList[listLen++] = this;
}

movable_Fl_Window::~movable_Fl_Window()
{
    for (int i = 0; i < listLen; i++) {
        if (windowList[i] == this) {
            windowList[i] = windowList[listLen - 1];
            listLen--;
            break;
        }
    }
}

int movable_Fl_Window::handle(int event)
{
    int ret;

#ifdef DEBUG
    cout << "movable_Fl_Window::handle()" << endl;
#endif
    if (event == FL_HIDE) {
        if (hideCallback != NULL)
            hideCallback();
    }
    else if (event == FL_SHOW) {
        if (showCallback != NULL)
            showCallback();
    }
    if (handleCallback != NULL)
        ret = handleCallback(event);
    else
        ret = 0;

    if (! ret)
        ret = Fl_Window::handle(event);
    return ret;
}

void movable_Fl_Window::resize(int x, int y, int w, int h)
{
#ifdef DEBUG
    cout << "movable_Fl_Window::resize(" << x << ", " << y << ", "
        << w << ", " << h << ')' <<endl;
#endif
    Fl_Window::resize(x, y, w, h);
    if (resizeCallback != NULL)
        (*resizeCallback)(x, y, w, h);  // Post-processing
}

void movable_Fl_Window::setShowCallback(void (* cb)())
{
    showCallback = cb;
}

void movable_Fl_Window::setHideCallback(void (* cb)())
{
    hideCallback = cb;
}

void movable_Fl_Window::setResizeCallback(void (* cb)(int, int, int, int))
{
    resizeCallback = cb;
}

void movable_Fl_Window::setHandleCallback(int (* cb)(int))
{
    handleCallback = cb;
}

