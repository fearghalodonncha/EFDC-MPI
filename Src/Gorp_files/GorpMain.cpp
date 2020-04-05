
#include "Gorp.h"

// Invoke with three arguments
// gorp cell.inp stem 25
int main(int argc, char **argv)
{
	Gorp g;
    cout << "Gorp begin " <<  "   "<<  endl;

    if (argc<4) {
        cerr << argv[0] << " CELLNAME.INP OutStem NActive" << endl;
        exit(-1);
    }

    g.doIt(argv[1], argv[2], atoi(argv[3]));  // min n

	return 0;
}
