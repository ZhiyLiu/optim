#include <stdio.h>
#include "M3DFigureTreeNode.h"

M3DFigureTreeNode::M3DFigureTreeNode(int _figureId)
{
    parent = NULL;
    sibling = NULL;
    figureId = _figureId;

    blendAmount = 0.0;
    blendExtent = 0.0;

	attachment = UNATTACHED;
    visited = false;
}

M3DFigureTreeNode::~M3DFigureTreeNode()
{
    int i;

    for(i = 0; i < children.size(); i++)
    {
        if(children[i] != NULL)
            delete children[i];
    }

    for(i = 0; i < links.size(); i++)
    {
        if(links[i] != NULL)
            delete links[i];
    }
}

M3DFigureTreeNode * M3DFigureTreeNode::copyPtr()
{
    M3DFigureTreeNode * retNode;
    int linkCount,
        childCount,
        i;
    M3DPrimitiveLinkInfo * oldLink,
                         * newLink;
    M3DFigureTreeNode * childNode;


    childCount = children.size();

    retNode = new M3DFigureTreeNode;

    retNode->figureId = figureId;
    retNode->blendAmount = blendAmount;
    retNode->blendExtent = blendExtent;
    retNode->attachment = attachment;

    // If the parent is not NULL (we aren't the root),
    // this will be set right in addChild()
    retNode->parent = NULL;
    retNode->sibling = NULL;

    linkCount = links.size();
    for(i = 0; i < linkCount; i++)
    {
        oldLink = links[i];
        if(oldLink == NULL)
            newLink = NULL;
        else
        {
            newLink = new M3DPrimitiveLinkInfo;
            newLink->primitiveId = oldLink->primitiveId;
            newLink->u = oldLink->u;
            newLink->v = oldLink->v;
            newLink->t = oldLink->t;
        }

        retNode->addLink(newLink);
    }

    for(i = 0; i < childCount; i++)
    {
        childNode = children[i];
        if(childNode != NULL)
            retNode->addChild(childNode->copyPtr());
    }

    return retNode;
}

void M3DFigureTreeNode::print(std::ostream & out, int depth)
{
    int childCount = children.size();
    int i;

	out << '\t';
    for(i = 0; i < depth; i++)
        out << "  ";

    out << figureId;
	out << " (";
	switch (getAttachmentMode())
	{
		case UNATTACHED:
			out << "main figure";
			break;
		case INDENT:
			out << "indentation";
			break;
		case PROTRUDE:
			out << "protrusion";
			break;
	}
	out << ")\n";

    for(i = 0; i < childCount; i++)
        if(children[i] != NULL)
            children[i]->print(out, depth + 1);
}

M3DFigureTreeNode * M3DFigureTreeNode::getChild(int id)
{
    if(id < 0 || id >= children.size())
        return NULL;

    return children[id];
}

M3DPrimitiveLinkInfo * M3DFigureTreeNode::getLink(int id)
{
    if(id < 0 || id >= links.size())
        return NULL;

    return links[id];
}

// Clears visited flag for this node and its subtree
void M3DFigureTreeNode::clearAllVisited()
{
    M3DFigureTreeNode * child;
    int childCount,
        i;


    visited = false;

    childCount = children.size();
    for(i = 0; i < childCount; i++)
    {
        child = children[i];
        if(child != NULL)
            child->clearAllVisited();
    }
}

// Counts the figures in this node and its subtree
int M3DFigureTreeNode::countFigures() const
{
    int figureCount;

    figureCount = 0;
	countTreeFigures(figureCount);
	return figureCount;
}

// Called recursively by M3DFigureTreeNode::countFigures()
void M3DFigureTreeNode::countTreeFigures(int & figureCount) const
{
    M3DFigureTreeNode * child;
    int childCount, i;

    figureCount++;

    childCount = children.size();
    for(i = 0; i < childCount; i++)
    {
        child = children[i];
        if(child != NULL)
            child->countTreeFigures(figureCount);
    }
}

// Decrements all figureId's greater than cutoffId in the tree
void M3DFigureTreeNode::decrement(int cutoffId)
{
    M3DFigureTreeNode * child;
    int childCount,
        i;


	if (figureId > cutoffId)
		figureId--;

    childCount = children.size();
    for(i = 0; i < childCount; i++)
    {
        child = children[i];
        if(child != NULL)
            child->decrement(cutoffId);
    }
}

void M3DFigureTreeNode::addChild(M3DFigureTreeNode * child)
{
    int childCount;
    M3DFigureTreeNode * prevChild;


    if(child == NULL)
        return;

    if(findNode(child->getFigureId()) != NULL)
    {
        printf("Error adding figure: figure is already in the tree.\n");
        return;
    }

    child->parent = this;
    child->sibling = NULL;

    childCount = children.size();
    if(childCount != 0)
    {
        prevChild = children[childCount - 1];
        if(prevChild != NULL)
            prevChild->sibling = child;
    }

    children.insert(children.end(), child);
}

// AGG: This function was not in use, so I commented it out
/*void M3DFigureTreeNode::removeChild(int childId)
{
    M3DFigureTreeNode * child;
    M3DFigureTreeNode * prevChild;
    int childCount;


    childCount = children.size();
    if(childId >= 0 && childId < childCount)
    {
        // First fix the sibling pointers for previous child
        if(childId > 0)
        {
            prevChild = children[childId - 1];
            if(childId < childCount - 2)
                prevChild->sibling = children[childId + 1];
            else
                prevChild->sibling = NULL;
        }

        // Now remove child from list
        child = children[childId];

        if(child != NULL)
            delete child;

        children.erase(children.begin() + childId);
    }
}*/

// Removes the child but doesn't delete it (returns the removed child)

M3DFigureTreeNode * M3DFigureTreeNode::popChild(int childId)
{
    M3DFigureTreeNode * child;
    M3DFigureTreeNode * prevChild;
    int childCount;


    childCount = children.size();
    if(childId >= 0 && childId < children.size())
    {
        // First fix the sibling pointers for previous child
        if(childId > 0)
        {
            prevChild = children[childId - 1];
            if(childId < childCount - 2)
                prevChild->sibling = children[childId + 1];
            else
                prevChild->sibling = NULL;
        }

        child = children[childId];

        children.erase(children.begin() + childId);

        return child;
    }

    return NULL;
}

void M3DFigureTreeNode::addLink(M3DPrimitiveLinkInfo * link)
{
    links.insert(links.end(), link);
}

void M3DFigureTreeNode::removeLink(int linkId)
{
    M3DPrimitiveLinkInfo * link;


    if(linkId >= 0 && linkId < links.size())
    {
        link = links[linkId];

        if(link != NULL)
            delete link;

        links.erase(links.begin() + linkId);
    }
}

// Finds the first node in the tree with figureId
M3DFigureTreeNode * M3DFigureTreeNode::findNode(int testFigureId)
{
    M3DFigureTreeNode * treeNode;
    int childCount;
    int i;


    if(figureId == testFigureId)
        return this;

    childCount = children.size();
    for(i = 0; i < childCount; i++)
    {
        treeNode = children[i]->findNode(testFigureId);
        if(treeNode != NULL)
            return treeNode;
    }

    return NULL;
}


