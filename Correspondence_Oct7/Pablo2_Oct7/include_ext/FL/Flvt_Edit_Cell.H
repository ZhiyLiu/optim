#ifndef FLVT_EDIT_CELL_H
#define FLVT_EDIT_CELL_H


/*  This class defines cells to be used by Flvt_Edit.
    Cell width can be explicitly set when a cell is created,
    or the global cell width will be used, unless it has
    not been defined, in which case the string length
    will be used.  When a cell is edited, it is the caller's
    responsibility to be certain that a longer string than
    the cell width is not inserted.  The cell width is one
    more than the possible number of characters that can be
    inserted.  An optional callback can be set with a unique
	tag; when the contents of the cell is changed, the callback
	function will be called with the tag as argument.
*/

class Flvt_Edit;

class Flvt_Edit_Cell {
	friend class Flvt_Edit;

	public:
		Flvt_Edit_Cell();
		Flvt_Edit_Cell(const char * value, int maxWidth = 0);
		virtual ~Flvt_Edit_Cell();

		void setCommonWidth(int maxWidth);
		int commonWidth() { return global_width; }
		int cellWidth() { return w; }

		void set(const char * value);
		void get(char * value);
		const char * textPtr() { return v; }

		void setCellChangeCallback(void (* fn)(void *, int), int tag);
		void setCellChangeObject(void * classObj) { ptr = classObj; }
		bool changed();

	private:
		char * v;
		int w;
		bool edited;

		static int global_width;

		void change(const char * value) {
			set(value);
			edited = true;
			if (callback)
				(* callback)(ptr, data);
		}
		void (* callback)(void *, int);
		int data;
		static void * ptr;	// Only one for all cells
};



/*  The following function is provided to create a pop-up, stand-alone
    Flvt_Edit widget.  Cells is an 1D array containing the 2D grid of
    editable cells with width of cols and height of rows.  Arguments
    x, y, w, and h and the usual FLTK window position and size parameters.
    The default width specifies the width of all columns.  If specified,
    col0Width overrides this width for the first column only.  The default
    width specifies the width of all columns.  If specified, col0Width
    overrides this width for the first column only.  Setting hscroll to
	true adds a horizontal scrollbar when needed.  The optional header
    specifies a single-line caption and label is the string to be placed
    in the window's banner.
*/

extern void flvt_popup_cell_editor(Flvt_Edit_Cell * cells, int rows, int cols,
					 int x, int y, int w, int h, int defaultWidth, int col0Width = 0,
					 bool hscroll = false, const char * header = NULL,
					 const char * label = NULL);


/*  The following functions embed an Flvt_Edit widget in a pre-existing
    window or group.  Cells is an 1D array containing the 2D grid of editable
    cells with width of cols and height of rows.  The first argument is the
    enclosing FLTK object.  The default width specifies the width of all columns.
    If specified, col0Width overrides this width for the first column only.
	Setting hscroll to true adds a horizontal scrollbar when needed.  The
    optional header specifies a single-line caption. The function returns a
    pointer to the Flvt_Edit object that is created.  If the callback function is
	specified, it will be called with the cell number as the second argument,
	every time a cell's value is changed.  The first argument will be the
	optionally specified classPtr, usually a "this" pointer.

    These functions add the Flvt_Edit object to the specified window or group.
    The caller will thus have to remove it from the window or group, before
    deleting the Flvt_Edit object.
*/

extern Flvt_Edit * flvt_cell_editor(class Fl_Window * win, Flvt_Edit_Cell * cells,
									int rows, int cols, int defaultWidth,
									int col0Width = 0, bool hscroll = false,
									const char * header = NULL,
									void (* changedCellCallback)(void *, int) = NULL,
									void * classPtr = NULL);

extern Flvt_Edit * flvt_cell_editor(class Fl_Group * grp, Flvt_Edit_Cell * cells,
									int rows, int cols, int defaultWidth,
									int col0Width = 0, bool hscroll = false,
									const char * header = NULL,
									void (* changedCellCallback)(void *, int) = NULL,
									void * classPtr = NULL);


#endif
