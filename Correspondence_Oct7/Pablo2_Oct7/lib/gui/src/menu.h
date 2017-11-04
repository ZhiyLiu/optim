// from PLUNC 6.4.0, based on modified fltk 1.5

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_Value_Slider.H>
//#include <FL/Fl_Double_Slider.H>	// PLUNC customization
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Scrollbar.H>
//#include <FL/Fl_Value_Scrollbar.H>	// PLUNC customization
#include <FL/Fl_Dial.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_Pixmap.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <FL/fl_draw.H>
#include <FL/fl_ask.H>
#include <stdlib.h>


void init_menu(const char *);
void add_menu_item(const char *, int value);
void add_menu_item(const char *, int value, int disp_flag);
int do_menu();

