// from PLUNC 6.4.0, based on modified fltk 1.5

#include "menu.h"

#ifndef NULL
#define NULL 0
#endif

static Fl_Menu_Button   *dummy;
Fl_Menu_Item            *menu_list;
static int              choice;

static void
dummy_cb(Fl_Widget *widget, void *ptr)
{   long val;
    
    val = (long)ptr;
    choice = val;
}

void init_menu(const char *label)
{   int ix, iy;
    ix = Fl::event_x_root();
    iy = Fl::event_y_root();
    dummy = new Fl_Menu_Button(ix, iy, 1, 1);
    menu_list = NULL;
    dummy->menu(menu_list);
    choice = -1;
    dummy->add(label, 0, dummy_cb, (void *)-1, FL_MENU_DIVIDER);
}

void add_menu_item(const char *label, int value, int disp_flag)
{   int flag = 0;

    if (disp_flag == -3) flag = 0;
//    if (disp_flag == -2) flag = FL_MENU_DOUBLE_DIVIDER;
    if (disp_flag == -2) flag = FL_MENU_DIVIDER;
    else if (disp_flag == -1) flag = FL_MENU_INACTIVE;
    else if (disp_flag == 0) flag = FL_MENU_TOGGLE;
    else if (disp_flag == 1) flag = FL_MENU_TOGGLE|FL_MENU_VALUE;
    if (flag == 0) dummy->add(label, 0, dummy_cb, (void *)value);
    else dummy->add(label, 0, dummy_cb, (void *)value, flag);
}

void add_menu_item(const char *label, int value)
{
    dummy->add(label, 0, dummy_cb, (void *)value);
}

int do_menu()
{
    dummy->popup();
    return(choice);
}
