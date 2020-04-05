#include "Gorp.h"


Gorp::~Gorp(void)
{
}


// Cell.inp file MUST have correct numbers at the top!!!  In that format...

//C Cell.inp file,  391 columns and  182 rows @           07/03/2008 18:08:32 
//C Project: Test1D
//C             1         2         3         4         5         6         7         8         9         0
//C    123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
//131  000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
void Gorp::loadTypes(char *fname) {
	FILE *f = fopen(fname, "r");
	if (!f)
		exit(-1);
	char buff[1024];
	char dummyname[1024];
	fgets(buff, 1024, f);
	sscanf(buff, "C %s file, %d columns and %d rows", dummyname, &mWidth, &mHeight);
	fgets(buff, 1024, f);
	fgets(buff, 1024, f);
	fgets(buff, 1024, f);
    int wtot = 0;
	mType = new int *[mHeight];
	for (int i=0; i<mHeight; i++)
		mType[i] = new int[mWidth];

    while (wtot < mWidth) {
        bool first = true;
        int wcurrent = 0;
        for (int j=mHeight-1; j>=0; j--) {
		    fgets(buff, 1024, f);
		    char *p = buff+5;
            int i;
		    for (i=0; *p >= '0' && *p <= '9'; i++, p++) {
			    mType[j][i+wtot] = *p-'0';
		    }
            if (first) {
                wcurrent = i;
                first = false;
            }
	    }
        wtot += wcurrent;
    }

#if 0
	for (int j=mHeight-1; j>=0; j--) {
		fgets(buff, 1024, f);
		char *p = buff+5;
		for (int i=0; i<mWidth-120; i++, p++) {
			mType[j][i+120] = *p-'0';
		}
	}
#endif

	fclose(f);

    ofstream cellout;
    // char buff[1024];
    sprintf(buff, "%s.txt", fname);
    cellout.open(buff);
    for (int j=mHeight-1; j>=0; j--) {
        for (int i=0; i<mWidth; i++) {
            cellout.put(char(mType[j][i])+'0');
        }
        cellout << endl;
    }
    cellout.close();
}
