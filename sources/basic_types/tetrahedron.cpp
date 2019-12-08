#include "tetrahedron.h"

void Tetrahedron::TE(int pos, vector<int>& e)
{
    e.assign(2,0);
    switch(pos)
    {
    case(0):
        e[0]=abs(vertices[0]); e[1]=abs(vertices[1]);
        break;
    case(1):
        e[0]=abs(vertices[0]); e[1]=abs(vertices[2]);
        break;
    case(2):
        e[0]=abs(vertices[0]); e[1]=abs(vertices[3]);
        break;
    case(3):
        e[0]=abs(vertices[1]); e[1]=abs(vertices[2]);
        break;
    case(4):
        e[0]=abs(vertices[1]); e[1]=abs(vertices[3]);
        break;
    case(5):
        e[0]=abs(vertices[2]); e[1]=abs(vertices[3]);
        break;
    default:
        cerr<<"[TE] wrong edge position. must be from 0 to 5."<<endl;
        int a; cin>>a;
        break;
    }
    sort(e.begin(),e.end());
}
