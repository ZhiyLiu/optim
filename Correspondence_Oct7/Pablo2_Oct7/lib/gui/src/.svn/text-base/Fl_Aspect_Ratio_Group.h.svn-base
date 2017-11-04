#ifndef FL_ASPECT_RATIO_GROUP
#define FL_ASPECT_RATIO_GROUP


class Fl_Aspect_Ratio_Group : public Fl_Group
{
public:
	Fl_Aspect_Ratio_Group(int x, int y, int w, int h, const char *label = 0) :
        Fl_Group(x, y, w, h, label) { }

    // Make sure window always has the same aspect ratio
    void resize(int x, int y, int w, int h);

private:

	static bool fullscreen;
	static bool last_status;
};


#endif
