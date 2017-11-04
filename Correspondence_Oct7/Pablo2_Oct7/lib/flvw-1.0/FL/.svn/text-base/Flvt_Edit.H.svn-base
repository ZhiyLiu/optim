#ifndef FLVT_EDIT_H
#define FLVT_EDIT_H


/*  This is a class for creating an editable table widget.  An
	Flvt_Edit object is intended to be embedded in one of the
	FLTK Window objects.  The containing window will typically
	be 20 pixels in X and 45 pixels in Y wider than the widget.
	The label field of the constructor can be used to specify a
	single, title cell across the top of the widget.  After
	construction, initialize() should be called to provide the
	array of cells to be displayed.  The left-most column cannot
	be edited, and is intended to be used for line numbers or
	names.
*/

class Flvt_Edit_Cell;

class Flvt_Edit : public Flv_Table
{
	public:
		Flvt_Edit(int X, int Y, int W, int H, const char * l = 0);
		virtual ~Flvt_Edit() {};
		void initialize(Flvt_Edit_Cell * cells, int rows) {
			bufs = cells;
			nrows = rows;
		}
		bool changed() { return changeFlag; }

		void end_edit(void);

	protected:
		int handle(int event);
		//	Required for editing
		virtual void save_editor( Fl_Widget *e, int R, int C );
		virtual void load_editor( Fl_Widget *e, int R, int C );
		virtual void position_editor( Fl_Widget *e, int x, int y,
			int w, int h, Flv_Style &s );

		//	Required for drawing
		virtual void draw_cell( int Offset, int &X, int &Y,
			int &W, int &H, int R, int C );

		const char * get_value(int R, int C);

		Flvt_Edit_Cell * bufs;
		int nrows;
		bool changeFlag;
};


#endif

