
#include <stdio.h>
#include <FL/Flv_Table.H>
#include <FL/Flve_Input.H>
#include <FL/Fl_Window.H>
#include <FL/fl_draw.H>
#include <FL/Enumerations.H>
#include <FL/Flvt_Edit.H>
#include <FL/Flvt_Edit_Cell.H>




Flvt_Edit * flvt_cell_editor(Fl_Window *win, Flvt_Edit_Cell * cells, int rows, int cols,
					int defaultWidth, int col0Width, bool hscroll, const char * header,
					void (* changedCellCallback)(void *, int), void * classPtr)

{
//	Flv_Style s;

	Flvt_Edit *e = new Flvt_Edit(10, 35, win->w() - 20, win->h() - 45, header );

	if (! hscroll)
		((Flv_List *) e)->has_scrollbar(FLVS_VERTICAL);

	e->end();
	Flve_Input *input = new Flve_Input( 0, 0, 0, 0 );
	win->end();

	input->hide();
	input->owner = e;
//	e->get_default_style(s);
	win->resizable(e);

//	e->edit_when(FLV_EDIT_ALWAYS);
	e->edit_when(FLV_EDIT_AUTOMATIC);
//	e->edit_when(FLV_EDIT_MANUAL);

	e->global_style.editor(input);

	e->select_locked(false);
	e->global_style.font_size(14);
	e->global_style.x_margin(5);
	e->global_style.locked(false);

	e->rows(rows);
	e->cols(cols);
	e->feature(FLVF_DIVIDERS|FLVF_PERSIST_SELECT);

	e->global_style.width(defaultWidth);
	e->global_style.height(15);
	e->global_style.font_size(12);
	input->textsize(12);

	e->col_style[0].align(FL_ALIGN_LEFT);
	if (col0Width > 0)
		e->col_style[0].width(col0Width);
	e->col_style[1].align(FL_ALIGN_LEFT);

	e->initialize(cells, rows);
	if (changedCellCallback)
		for (int i = 0; i < rows*cols; i++)
			cells[i].setCellChangeCallback(changedCellCallback, i);
	cells[0].setCellChangeObject(classPtr);

	e->row_style[-1].align(FL_ALIGN_CLIP);

	e->col_style[-1].align(FL_ALIGN_CLIP);

	e->col_style[0].locked(true);

	return e;
}


Flvt_Edit * flvt_cell_editor(Fl_Group *grp, Flvt_Edit_Cell * cells, int rows, int cols,
					int defaultWidth, int col0Width, bool hscroll, const char * header,
					void (* changedCellCallback)(void *, int), void * classPtr)
{
//	Flv_Style s;

	Flvt_Edit *e = new Flvt_Edit(grp->x() + 5, grp->y() + 5, grp->w() - 10, grp->h() - 10, header );

	if (! hscroll)
		((Flv_List *) e)->has_scrollbar(FLVS_VERTICAL);

	e->end();
	Flve_Input *input = new Flve_Input( 0, 0, 0, 0 );
	grp->add(e);	// AGG: Testing
	grp->end();

	input->hide();
	input->owner = e;
//	e->get_default_style(s);
	grp->resizable(e);

//	e->edit_when(FLV_EDIT_ALWAYS);
	e->edit_when(FLV_EDIT_AUTOMATIC);
//	e->edit_when(FLV_EDIT_MANUAL);

	e->global_style.editor(input);

	e->select_locked(false);
	e->global_style.font_size(14);
	e->global_style.x_margin(5);
	e->global_style.locked(false);

	e->rows(rows);
	e->cols(cols);
	e->feature(FLVF_DIVIDERS|FLVF_PERSIST_SELECT);

	e->global_style.width(defaultWidth);
	e->global_style.height(15);
	e->global_style.font_size(12);
	input->textsize(12);

	e->col_style[0].align(FL_ALIGN_LEFT);
	if (col0Width > 0)
		e->col_style[0].width(col0Width);
	e->col_style[1].align(FL_ALIGN_LEFT);

	e->initialize(cells, rows);
	if (changedCellCallback)
		for (int i = 0; i < rows*cols; i++)
			cells[i].setCellChangeCallback(changedCellCallback, i);
	cells[0].setCellChangeObject(classPtr);

	e->row_style[-1].align(FL_ALIGN_CLIP);

	e->col_style[-1].align(FL_ALIGN_CLIP);

	e->col_style[0].locked(true);
#ifdef DEBUG
	printf("e = 0x%x\n", e);
#endif
	return e;
}


void flvt_popup_cell_editor(Flvt_Edit_Cell * cells, int rows, int cols,
					 int x, int y, int w, int h, int defaultWidth, int col0Width,
					 bool hscroll, const char *header, const char * label) 
{
	Fl_Window *win = new Fl_Window(x, y, w, h, label);

	Flvt_Edit * e = flvt_cell_editor(win, cells, rows, cols,
		defaultWidth, col0Width, hscroll, header);

	win->show(0, NULL);
}

