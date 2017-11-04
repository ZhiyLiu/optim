#ifndef M3D_FIGURE_TREE_NODE_H
#define M3D_FIGURE_TREE_NODE_H

#include <vector>
#include <iostream>


struct M3DPrimitiveLinkInfo
{
    double u, v, t;		// Figural coord's of linked atom's attachment to main figure
    int primitiveId;	// ID of linked atom in attached subfigure
};


/*  See P3DControl::attachSubfigure() for an example of the use of this class.
    Note that a node with a NULL parent is the root node of a tree.  By default
	the blendAmount and blendExtent in a root node are set to 0.0.  Links are
	placed in the child figure that is linking to a parent.
*/

class M3DFigureTreeNode
{

public:

	enum SubfigureAttachment_t
	{
		// These are assigned values to simplify the user interface.
		// See P3DUserInterfaceCallback::setAttachSubfigureMode(int mode).
		UNATTACHED = 0,
		PROTRUDE = 1,
		INDENT = 2
	};

    M3DFigureTreeNode(int _figureId = -1);
    ~M3DFigureTreeNode();

    // Copies the entire (sub)tree rooted at this node,
    // and returns a pointer to the copied root
    M3DFigureTreeNode * copyPtr();

    void print(std::ostream & out = std::cout, int depth = 0);

    int getChildCount() const { return children.size(); }
    int getLinkCount() const { return links.size(); }
    int countFigures() const;

    int getFigureId() const { return figureId; }
    double getBlendAmount() const { return blendAmount; }
    double getBlendExtent() const { return blendExtent; }
    M3DFigureTreeNode * getParent() { return parent; }
    M3DFigureTreeNode * getChild(int id);
    M3DFigureTreeNode * getSibling() { return sibling; }
    M3DPrimitiveLinkInfo * getLink(int id);
    bool getVisited() const { return visited; }
    SubfigureAttachment_t getAttachmentMode() const { return attachment; }

    void setFigureId(int _figureId) { figureId = _figureId; }
    void setBlendAmount(double _blendAmount) { blendAmount = _blendAmount; }
    void setBlendExtent(double _blendExtent) { blendExtent = _blendExtent; }
    void setParent(M3DFigureTreeNode * _parent) { parent = _parent; }
    void setVisited(bool _visited) { visited = _visited; }
    void setAttachmentMode(SubfigureAttachment_t mode) { attachment = mode; }

    // Clears visited flag for this node and its subtree
    void clearAllVisited();

    void addChild(M3DFigureTreeNode * child);

    // Deletes the subtree rooted at the childId
//   void removeChild(int childId);

    // Removes the child but doesn't delete it (returns the removed child)
    M3DFigureTreeNode * popChild(int childId);

    void addLink(M3DPrimitiveLinkInfo * link);
    void removeLink(int linkId);

	// Function used when figures are deleted
	void decrement(int cutoffId);

    // Finds the first node in the tree with figureId
    M3DFigureTreeNode * findNode(int testFigureId);

private:

    int figureId;

    // Parameters defining how the figure and parent surfaces are blended
    double blendAmount;
    double blendExtent;

    // Pointer up to parent node
    M3DFigureTreeNode * parent;

    // Pointer across to next sibling
    M3DFigureTreeNode * sibling;

    // List of pointers to children nodes
    std::vector<M3DFigureTreeNode *> children;

    // Information for primitives that are linked to the parent figure
    std::vector<M3DPrimitiveLinkInfo *> links;

    // Flag for traversal
    bool visited;

	// Type of attachment
	SubfigureAttachment_t attachment;

    void countTreeFigures(int & figureCount) const;
};

#endif

