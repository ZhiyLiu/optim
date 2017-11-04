
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <iostream>

#define LINE_WIDTH	2
#define DOT_WIDTH	1



/*   ^
     | \
     |  \
	 +---+
*/
void drawTriangle_1(int x, int y)
{
	fl_line(x + 3, y + 3, x + 12, y + 18);	// hypotenuse
	fl_line(x + 12, y + 18, x + 3, y + 18);	// base
	fl_line(x + 3, y + 18, x + 3, y + 3);	// vertical
}

/*       ^
       / |
      /  |
	 +---+
*/
void drawTriangle_2(int x, int y)
{
	fl_line(x + 19, y + 3, x + 9, y + 18);	// hypotenuse
	fl_line(x + 9, y + 18, x + 19, y + 18);	// base
	fl_line(x + 19, y + 18, x + 19, y + 3);	// vertical
}

/*   +---->
     |  /
     +/
*/
void drawTriangle_3(int x, int y, bool longer = false)
{
	int x0;

	if (longer)
		x0 = 4;
	else
		x0 = 5;
	fl_line(x + x0, y + 3, x + 19, y + 3);	// vertical
	fl_line(x + 19, y + 3, x + x0, y + 12);	// hypotenuse
	fl_line(x + x0, y + 12, x + x0, y + 3);	// base
}

/*   <----+
       \  |
         \+
*/
void drawTriangle_4(int x, int y)
{
	fl_line(x + 17, y + 3, x + 3, y + 3);	// vertical
	fl_line(x + 3, y + 3, x + 17, y + 12);	// hypotenuse
	fl_line(x + 17, y + 12, x + 18, y + 3);	// base
}

/*   +\
     |  \
     +---->
*/
void drawTriangle_5(int x, int y)
{
	fl_line(x + 5, y + 18, x + 19, y + 18);	// vertical
	fl_line(x + 19, y + 18, x + 5, y + 9);	// hypotenuse
	fl_line(x + 5, y + 9, x + 5, y + 18);	// base
}

void drawRotateZLeft90(Fl_Label *label, int x, int y, int w, int h, Fl_Align align)
{
	fl_color(4);
	fl_line_style(FL_SOLID, LINE_WIDTH);
	drawTriangle_1(x, y);

	fl_line_style(FL_DOT, DOT_WIDTH, NULL);
	drawTriangle_3(x, y);
	fl_line_style(FL_SOLID);
}

void drawRotateZRight90(Fl_Label *label, int x, int y, int w, int h, Fl_Align align)
{
	fl_color(4);
	fl_line_style(FL_SOLID, LINE_WIDTH);
	drawTriangle_2(x, y);

	fl_line_style(FL_DOT, DOT_WIDTH, NULL);
	drawTriangle_4(x - 1, y);
	fl_line_style(FL_SOLID);
}

void drawRotateX180(Fl_Label *label, int x, int y, int w, int h, Fl_Align align)
{
	fl_color(3);
	fl_line_style(FL_SOLID, LINE_WIDTH);
	drawTriangle_1(x, y);

	fl_line_style(FL_DOT, DOT_WIDTH, NULL);
	drawTriangle_2(x, y);
	fl_line_style(FL_SOLID);
}

void drawRotateY180(Fl_Label *label, int x, int y, int w, int h, Fl_Align align)
{
	fl_color(5);
	fl_line_style(FL_SOLID, LINE_WIDTH);
	drawTriangle_5(x, y);

	fl_line_style(FL_DOT, DOT_WIDTH, NULL);
	drawTriangle_3(x, y, true);
	fl_line_style(FL_SOLID);
}

// A dispatcher is required because there can only be one "free label type".
// Because there is no other easy way to distinguish the four buttons, the
// alignment code is used.  This means, it cannot be changed in P3DUserInterface.fl.
void drawRotateLabel(Fl_Label *label, int x, int y, int w, int h, Fl_Align align)
{
	switch (align) {
		case 16:	drawRotateZLeft90(label, x + 2, y, w, h, align);
					break;
		case 20:	drawRotateX180(label, x + 2, y, w, h, align);
					break;
		case 18:	drawRotateY180(label, x + 2, y, w, h, align);
					break;
		case 24:	drawRotateZRight90(label, x + 2, y, w, h, align);
					break;
		default:	std::cout << "Invalid call to drawRotateLabel()\n";
					break;
	}
}

