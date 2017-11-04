
// This class was derived from an earlier test program and included in
// the Flvt library by A. G. Gash, 23 July 2001.  The original comments
// are given below.

//	======================================================================
//	File:    edit.cxx
//	Version: 0.1.0
//	Started: 01/03/00
//
//	Copyright (C) 1999 Laurence Charlton
//
//	Description:
//	Integrate editing.  DO NOT USE THIS FILE!!!! IT IS ONLY FOR
//	INTERNAL USE.  FINAL EDITING MOST LIKELY WILL BEAR NO SIMILARITIES
//	WITH THIS CODE.  I AM USING IT TO SEE WHAT I CAN DO...
//
//	Notes:
//	** Fix callback_when to work
//	** Fix edit_when to work
//	** Fix locked to work
//	** Fix select_locked to work
//	** Test out of cell editors
//	Test large editors(drop-down)
//	Test other editors
//	Fix Flv_List to do editing
//	** Check if we need to set parent for editor
//	======================================================================

#include <FL/Flv_Table.H>
#include <FL/Flve_Input.H>
#include <FL/Fl_Window.H>
#include <FL/fl_draw.H>
#include <FL/Enumerations.H>
#include <FL/Flve_Check_Button.H>
#include <stdio.h>
#include <string.h>
#include <FL/Flvt_Edit.H>
#include <FL/Flvt_Edit_Cell.H>



Flvt_Edit::Flvt_Edit(int X, int Y, int W, int H, const char * l) : Flv_Table(X, Y, W, H, l)
{
#ifdef DEBUG
	printf("Flvt_Edit::Flvt_Edit()\n");
#endif
	changeFlag = false;
}

/*Flvt_Edit(const Flvt_Edit & edit)
{
	changeFlag = edit.changeFlag;
	nrows = edit.nrows;
	bufs = edit.bufs;
}*/

void Flvt_Edit::end_edit(void)
{
	Flv_Table::end_edit();
}

const char * Flvt_Edit::get_value(int R, int C)
{
	static char buf[40];
	*buf = '\0';
	if (R==-1 && C>-1)				//	Row header, A, B, C...
	{
		sprintf( buf, "  %c%c  ", C/26 + 'A'-1, (C%26)+'A' );
		if (*(buf+2)<'A')
			memmove( buf+2, buf+3, 4 );
	} else if (C==-1 && R>-1)	//	Column header 1, 2, 3...
		sprintf( buf, "%d", R );
	else if (R>-1 && C>-1) {	//	Normal cell from bufs
		int n = R + C*nrows;
		return bufs[n].textPtr();
	}
	return buf;
}

//	Required for manual editing
int Flvt_Edit::handle(int event)
{
#ifdef DEBUG
	printf("E:%d\n", event);
#endif
	fflush(stdout);
	if (event==FL_KEYBOARD)
	{
		changeFlag = true;
		if (Fl::event_key()==FL_F+2)	// Function key F2
		{
			start_edit();
			return 1;
		}
	}
	return Flv_Table::handle(event);
}

//	Required for editing
void Flvt_Edit::save_editor( Fl_Widget *e, int R, int C )	// Called when editing completes
{
	if (C > 0) {	// The left-most column cannot be edited
		int n = R + C*nrows;
		bufs[n].change(((Flve_Input *)e)->value());
	}
}

void Flvt_Edit::load_editor( Fl_Widget *e, int R, int C )	// Called when editing begins
{
	if (C > 0) {	// The left-most column cannot be edited
		int n = R + C*nrows;
		((Flve_Input *)e)->value(bufs[n].textPtr());
		((Flve_Input *)e)->position(((Flve_Input *)e)->size(), 0);
	}
}

//	Required for drawing
void Flvt_Edit::draw_cell( int Offset, int &X, int &Y, int &W, int &H, int R, int C )
{
	Flv_Style s;

#ifdef DEBUG
	printf("Flvt_Edit::draw_cell()\n");
#endif
	get_style(s, R, C);
	Flv_Table::draw_cell(Offset,X,Y,W,H,R,C);
	fl_draw(get_value(R,C), X-Offset, Y, W, H, s.align());
}

void Flvt_Edit::position_editor( Fl_Widget *e, int x, int y, int w, int h, Flv_Style &s )
{

//	Out of cell
//	e->resize( 10, 10, 200, 20 );

//	In cell
//	Flv_Table::position_editor(e,x+s.x_margin(),y,w-s.x_margin(),h,s);
	Flv_Table::position_editor(e,x,y,w,h,s);
}

