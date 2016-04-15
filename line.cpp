#include "line.h"
//=========================================================================
//
//                               Line  class
//
//=========================================================================
LINE::LINE ()
{
    Ncell = 0;
    Length = 0;
}

void LINE::Update ()
{
    double sPointer = 0.;
    for (unsigned i = 0; i < Cell.size (); i++)
    {
        Cell[i].S = sPointer;
        sPointer = Cell[i].L + sPointer;
    }
    Ncell = Cell.size ();
    Length = sPointer;
}

void LINE::Append (ELEMENT elem)
{
    Cell.push_back (elem);
    elem.S=Length;
    Length += elem.L;
    Ncell+=1;
}
