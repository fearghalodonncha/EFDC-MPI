#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

using namespace std;

class Cut
{
public:
	int mPos;
	bool mHoriz;
};

class Layout
{
public:
	double mMaxCost;
	double *mCosts;
	int *mHCuts;
	int *mVCuts;
	int mNHCuts;
	int mNVCuts;
	int mNPartitions;
	int mNActivePartitions;

	Layout() {
		mMaxCost = 0;
		mCosts = NULL;
		mHCuts = NULL;
		mVCuts = NULL;
		mNHCuts = 0;
		mNVCuts = 0;
		mNPartitions = 0;
		mNActivePartitions = 0;
	}

	void init(int nh, int nv) {
		mMaxCost = 0;
		mNHCuts = nh;
		mNVCuts = nv;
		mNPartitions = nh*nv;
		if (mCosts)
			delete [] mCosts;
		mCosts = new double[mNPartitions];
		if (mHCuts)
			delete [] mHCuts;
		mHCuts = new int[mNHCuts];
		if (mVCuts)
			delete [] mVCuts;
		mVCuts = new int[mNVCuts];
	}

    void moveHCut(int i, int d) {
        mHCuts[i]+=d;
        mHCuts[i+1]-=d;
    }

    void moveVCut(int i, int d) {
        mVCuts[i]+=d;
        mVCuts[i+1]-=d;
    }

	void init(int nh, int nv, int *hc, int *vc, double cost, int nactive) {
		init(nh, nv);
		memcpy(mHCuts, hc, nh*sizeof(int));
		memcpy(mVCuts, vc, nv*sizeof(int));
		mMaxCost = cost;
		mNActivePartitions = nactive;
	}

	void initVals(int target, int mint) {
		mMaxCost = 0;
	}

	Layout(const Layout& other) {
		mMaxCost = other.mMaxCost;
		mNPartitions = other.mNPartitions;
		mNActivePartitions = other.mNActivePartitions;
		mNHCuts = other.mNHCuts;
		mNVCuts = other.mNVCuts;
		mCosts = new double[mNPartitions];
		memcpy(mCosts, other.mCosts, mNPartitions*sizeof(double));
		mHCuts = new int[mNHCuts];
		memcpy(mHCuts, other.mHCuts, mNHCuts*sizeof(int));
		mVCuts = new int[mNVCuts];
		memcpy(mVCuts, other.mVCuts, mNVCuts*sizeof(int));
	}

    ~Layout() {
        if (mCosts)
            delete [] mCosts;
        if (mHCuts)
            delete [] mHCuts;
        if (mVCuts)
            delete [] mVCuts;
    }

	Layout& operator=(const Layout& other) {
		mMaxCost = other.mMaxCost;
		mNPartitions = other.mNPartitions;
		mNActivePartitions = other.mNActivePartitions;
		mNHCuts = other.mNHCuts;
		mNVCuts = other.mNVCuts;
		if (mCosts)
			delete [] mCosts;
		mCosts = new double[mNPartitions];
		memcpy(mCosts, other.mCosts, mNPartitions*sizeof(double));
		if (mHCuts)
			delete [] mHCuts;
		mHCuts = new int[mNHCuts];
		memcpy(mHCuts, other.mHCuts, mNHCuts*sizeof(int));
		if (mVCuts)
			delete [] mVCuts;
		mVCuts = new int[mNVCuts];
		memcpy(mVCuts, other.mVCuts, mNVCuts*sizeof(int));
		return *this;
	}

	inline bool operator > (const Layout& rhs) const {
		if (mMaxCost > rhs.mMaxCost)
			return true;
		if (mMaxCost < rhs.mMaxCost)
			return false;
		// if cost is equal...
		if (mNPartitions > rhs.mNPartitions)
			return false;
		if (mNPartitions < rhs.mNPartitions)
			return true;
		return false;
	}
};

class Partition
{
public:
	// X=2, L=3 means 2, 3, 4  So next partition would be 5, 3
	int mX, mY, mL, mW;
	double mCost;
    double mAreaWithGhostZone;
    double mGhostZone;
    double mNearLand;
    double mWater;
    double mLand;
	int mID;
	void set(int x, int y, int l, int w, int id) {
		mX = x;
		mY = y;
		mL = l;
		mW = w;
		mID = id;
	}

	inline bool operator == (const Partition& rhs) const {
		if (mX == rhs.mX && mY == rhs.mY && mL == rhs.mL && mW == rhs.mW)
			return true;
		return false;
	}

	inline Partition& operator = (const Partition& rhs) {
		mX = rhs.mX;
		mY = rhs.mY;
		mL = rhs.mL;
		mW = rhs.mW;
		mID = rhs.mID;
		mCost = rhs.mCost;
        mAreaWithGhostZone = rhs.mAreaWithGhostZone;
        mGhostZone = rhs.mGhostZone;
        mNearLand = rhs.mNearLand;
        mWater = rhs.mWater;
		return *this;
	}
};

class Gorp
{
public:
	// Gorp(void);
	~Gorp(void);

	int mWidth;
	int mHeight;

	enum CELLTYPE {LAND=0, NEARLAND=9, WATER=5};
	// enum {MINTHICK=25};
    int MINTHICKX;
    int MINTHICKY;

	int mTargetNActiveParts;
	int mNPartitions;
	int **mType;


	Partition *mPartitions;
	Partition *mPrevPartitions;

    int *mCutsHoriz;
    int *mCutsVert;

	int mNHoriz;
	int mNVert;
    int mTotalLand;
    int mTotalWater;

    double mTotalAreaWeight;
    double mTileAreaWeight;
    double mTileLandWeight;
    double mTileGhostWeight;
    double mTileWaterWeight;
    double mTileNearLandWeight;

    ofstream mLogOut;

	Gorp() {
		mWidth = 0;
		mHeight = 0;
		mType = NULL;
        mCutsHoriz = NULL;
        mCutsVert = NULL;
		mPartitions = NULL;
        mPrevPartitions = NULL;
		mNPartitions = 0;
		mTargetNActiveParts = 0;
        MINTHICKX = 0;
        MINTHICKY = 0;

        mTotalAreaWeight=0.00/203/131;  // multiplies total area;  approx 5ms times 2
        mTileAreaWeight=0.00/44/25;
        mTileGhostWeight=0.001/(4*44+4*25);
        mTileWaterWeight=0.007/44/25; 
        mTileLandWeight=0.00/44/25;
        mTileNearLandWeight=0.003/44/25;
    }

	void loadTypes(char *fname);

    // estimate the cost of a single tile
	inline double setCost(Partition &p) {
		p.mCost = 0;
        if (p.mY<0 || p.mY+p.mW > mHeight || p.mX<0 || p.mX+p.mL > mWidth || p.mL <= 0 || p.mW <= 0) {
            p.mCost = 1E50;
            return p.mCost;
        }
        double watersum=0;
        double nearlandsum=0;
        double landsum=0;
		for (int j=p.mY; j<p.mY+p.mW; j++) {
			for (int i=p.mX; i<p.mX+p.mL; i++) {
				if (mType[j][i] == WATER)
					watersum += 1.0;
				else if (mType[j][i] == NEARLAND)
					nearlandsum += 1;
                else if (mType[j][i] == LAND)
                    landsum += 1;
			}
		}
        p.mWater = watersum;
        p.mNearLand = nearlandsum;
        p.mGhostZone = 2*p.mW+2*p.mL;
        p.mAreaWithGhostZone = (p.mW+4)*(p.mL+4);
        p.mLand = landsum;

        double esttotwater = watersum/(p.mW*p.mL)*p.mAreaWithGhostZone;  // estimate water assuming ghost zone is similar

        // if no water, cost is 0
        if (watersum == 0 && nearlandsum == 0) {
            p.mCost = 0;
            return 0;
        }
		p.mCost = mTileGhostWeight*(2*p.mW+2*p.mL)+mTileWaterWeight*esttotwater+mTileNearLandWeight*nearlandsum+mTileLandWeight*landsum;  // local cost does not include area
		return p.mCost;
	}

	void setPartitionCount(int nh, int nv) {
		mNHoriz = nh;
		mNVert = nv;
		mNPartitions = (nh+1)*(nv+1);
		if (mPartitions)
			delete [] mPartitions;
		mPartitions = new Partition[mNPartitions];
	}

	// set all partitions based on cuts
    // does not set costs
	void setPartitions() {
		int npart = 0;
		for (int j=0; j<mNVert; j++) {
			int y0 = mCutsVert[j];
			int W = mCutsVert[j+1]-y0;
			for (int i=0; i<mNHoriz; i++) {
				int x0 = mCutsHoriz[i];
				int L = mCutsHoriz[i+1]-x0;
				mPartitions[npart].set(x0, y0, L, W, npart);
				npart++;
			}
		}	
	}

    // finds cost of each tile 
	void findStats(int &nactive, double &maxcost) {
		maxcost=0;
		nactive = 0;
        double activearea = 0;
		for (int i=0; i<mNPartitions; i++) {
			double c = setCost(mPartitions[i]);
			if (c>0) {
				nactive++;
				if (maxcost < c)
					maxcost = c;
                activearea += mPartitions[i].mL*mPartitions[i].mW;
			}
		}
        // now we add the area cost to max cost
        maxcost += activearea*mTotalAreaWeight;
    }

    // costs already calculated
    void findQuickPrevStats(int &nactive, double &maxcost) {
        maxcost = 0;
        nactive = 0;
        double activearea = 0;
		for (int i=0; i<mNPartitions; i++) {
			double c = mPrevPartitions[i].mCost;
			if (c>0) {
				nactive++;
				if (maxcost < c)
					maxcost = c;
                activearea += mPrevPartitions[i].mL*mPrevPartitions[i].mW;  //tile area 
			}
		}
        maxcost += activearea*mTotalAreaWeight;   //tile area * tile weight
    }

	// set all partitions using 0-based widths
	// counts are number of partitions
	void setPartitions(int nh, int nv, int *hc, int *hv) {
        if (nh != mNHoriz || !mCutsHoriz) {
            mNHoriz = nh;
            if (mNHoriz)
                delete [] mCutsHoriz;
            mCutsHoriz = new int[mNHoriz];
        }
        if (nv != mNVert || !mCutsVert) {
            mNVert = nv;
            if (mCutsVert)
                delete [] mCutsVert;
            mCutsVert = new int[mNVert];
        }
        for (int i=0; i<nh; i++)
            mCutsHoriz[i] = hc[i];
        for (int i=0; i<nv; i++)
            mCutsVert[i] = hv[i];
		if (mNPartitions != nh*nv) {
			mNPartitions = nh*nv;
			if (mPartitions)
				delete [] mPartitions;
			if (mPrevPartitions)
				delete [] mPrevPartitions;
			mPartitions = new Partition[mNPartitions];
			mPrevPartitions = new Partition[mNPartitions];
			int npart = 0;
			int y0=0;
			for (int j=0; j<mNVert; j++) {
				int W = hv[j]+MINTHICKY;
				int x0=0;
				for (int i=0; i<mNHoriz; i++) {
					int L = hc[i]+MINTHICKX;
					mPartitions[npart].set(x0, y0, L, W, npart);
					setCost(mPartitions[npart]);
					mPrevPartitions[npart] = mPartitions[npart];
					npart++;
					x0 += L;
				}
				y0 += W;
			}
		} else {
			int npart = 0;
			int y0=0;
			for (int j=0; j<mNVert; j++) {
				int W = hv[j]+MINTHICKY;
				int x0=0;
				for (int i=0; i<mNHoriz; i++) {
					int L = hc[i]+MINTHICKX;
					mPartitions[npart].set(x0, y0, L, W, npart);
					if (mPartitions[npart] == mPrevPartitions[npart])
						mPartitions[npart] = mPrevPartitions[npart];
					else {
						setCost(mPartitions[npart]);
						mPrevPartitions[npart] = mPartitions[npart];
					}
					npart++;
					x0 += L;
				}
				y0 += W;
			}
		}
	}

	void initPerm(int *v, int k, int n) {
		for (int i=0; i<k; i++)
			v[i] = 0;
		v[0] = n;
	}

	void displayPerm(int *v, int k) {
		for (int i=0; i<k; i++)
			cout << v[i] << " ";
		cout << endl;
	}

	// all need to sum to n
	// k is number of bins
	// start with n 0 0 0
	// n-1 1 0 0
	// 0 n 00
	// n-1 0 1
	// p points to highest activ index
	int nextPerm(int *v, int k, int n) {
		// always init with 0 as full as possible and empty into next one
		if (v[0]>0) {
			v[0]--;
			v[1]++;
			return 1;
		}
		// 0 bin is empty so decrement next non empty one and bump adjacent
		int j = 1;
		while (v[j] == 0)
			j++;
		if (j>=k-1)
			return 0;
		// found first nonzero column, so decrement and bump next one
		v[j+1]++;
		v[0] = v[j]-1;
		for (int i=1; i<=j; i++)
			v[i] = 0;
		return 1;
	}

        void findTotalLandWater() {
            mTotalLand = 0;
            mTotalWater = 0;
            for (int j=0; j<mHeight; j++)
                for (int i=0; i<mWidth; i++) {
                    if (mType[j][i] == LAND)
                        mTotalLand++;
                    else if (mType[j][i] == WATER)
                        mTotalWater++;
                }
                    

        }

        // calls findStats and calcs cost of each tile
        void findCost(Layout &l) {
            setPartitions(l.mNHCuts, l.mNVCuts, l.mHCuts, l.mVCuts);
            findStats(l.mNActivePartitions, l.mMaxCost);
            for (int i=0; i<l.mNPartitions; i++)
                l.mCosts[i] = mPartitions[i].mCost;
        }

        void findApproxLayout(int nh, int nv, Layout &l) {
	    int *hc = new int[nh];
	    int *vc = new int[nv];
            findTotalLandWater();
            //int wtarget = mTotalWater/nv;
            int nrow = 0;
            int wsum = 0;
            int prevrow = 0;
            int currentsum = 0;
            for (int i=0; i<nv-1; i++) {
                int wtarget = (i+1)*float(mTotalWater)/nv;
                while(currentsum<wtarget && nrow < mHeight) {
                    for (int j=0; j<mWidth; j++) {
                        if (mType[nrow][j] == WATER) {
                            currentsum++;
                        }
                    }
                    nrow++;
                }
                vc[i] = nrow-prevrow;
                wsum += vc[i];
                prevrow = nrow;
            }
            vc[nv-1] = mHeight-wsum;

            //wtarget = mTotalWater/nh;
            int ncol = 0;
            int hsum = 0;
            int prevcol = 0;
            currentsum = 0;
            for (int i=0; i<nh-1; i++) {
                int wtarget = (i+1)*float(mTotalWater)/nh;
                while(currentsum<wtarget) {
                    for (int j=0; j<mHeight; j++)
                        if (mType[j][ncol] == WATER)
                            currentsum++;
                    ncol++;
                }
                hc[i] = ncol-prevcol;
                hsum += hc[i];
                prevcol = ncol;
            }
            hc[nh-1] = mWidth-hsum;

            l.init(nh, nv, hc, vc, 0, 0);
            delete [] hc;
            delete [] vc;
        }

		double moveHCost(int h, int d, int maxactive) {
			for (int i=0; i<mNPartitions; i++)
				mPrevPartitions[i] = mPartitions[i];  // just set layout and costs
			for (int j=0; j<mNVert; j++) {  // go down vertically through column and squeeze horizontally
				int n = j*mNHoriz+h;        // index of partition in current row
				mPrevPartitions[n].mL += d;
                mPrevPartitions[n+1].mX += d;
                mPrevPartitions[n+1].mL -= d;
                setCost(mPrevPartitions[n]);
                setCost(mPrevPartitions[n+1]);
			}
            int nactive;
            double cost;
            findQuickPrevStats(nactive, cost);
            if (nactive!=maxactive)
                cost += 1E10*abs(1+nactive-maxactive);
			return cost;
		}

		double moveVCost(int v, int d, int maxactive) {
			for (int i=0; i<mNPartitions; i++)
				mPrevPartitions[i] = mPartitions[i];  // just set layout and costs
			for (int j=0; j<mNHoriz; j++) {  // go across horizontally through column and squeeze vertically
				int n = v*mNHoriz+j;         // index of partition in current column
				mPrevPartitions[n].mW += d;
                mPrevPartitions[n+mNHoriz].mY += d;
                mPrevPartitions[n+mNHoriz].mW -= d;
                setCost(mPrevPartitions[n]);
                setCost(mPrevPartitions[n+mNHoriz]);
			}
            int nactive;
            double cost;
            findQuickPrevStats(nactive, cost);
            if (nactive!=maxactive)
                cost += 1E10*abs(1+nactive-maxactive);
			return cost;
		}

		// need partitions set matching orig layout
		// use prevpartition to play games
		// this moves a single divider either + or - d in h or v direction - whichever SINGLE move has biggest gain
		// this changes lorig if a single move helps, and it does the best single move
        bool improvePartition(int d, Layout &lorig, int maxactive) {
			Layout l;
			bool hit = false;
			Layout best;
			double bestcost;

            findCost(lorig);
            if (lorig.mNActivePartitions != maxactive)
                lorig.mMaxCost += 1E10*abs(1+lorig.mNActivePartitions-maxactive);


            bestcost = lorig.mMaxCost;

            if (bestcost < 1E9) {
                mLogOut << "step " << d << " " << lorig.mNHCuts << " " << lorig.mNVCuts << " " << bestcost << endl;
            }

			for (int i=0; i<lorig.mNHCuts-1; i++) {
				l = lorig;
				double mcost;
				if ((mcost = moveHCost(i, d, maxactive)) < bestcost) {
					l.moveHCut(i, d);
					best = l;
					bestcost = mcost;
					hit = true;
                    if (bestcost < 1E10)
                        mLogOut << lorig.mNHCuts << " " << lorig.mNVCuts << " " << bestcost << endl;
				}
                l = lorig;
				if ((mcost = moveHCost(i, -d, maxactive)) < bestcost) {
					l.moveHCut(i, -d);
					best = l;
					bestcost = mcost;
					hit = true;
                    if (bestcost < 1E10)
                        mLogOut << lorig.mNHCuts << " " << lorig.mNVCuts << " " << bestcost << endl;
				}
			}

			for (int i=0; i<lorig.mNVCuts-1; i++) {
                l = lorig;
				double mcost;
				if ((mcost = moveVCost(i, d, maxactive)) < bestcost) {
					l.moveVCut(i, d);
					best = l;
					bestcost = mcost;
					hit = true;
                    if (bestcost < 1E10)
                        mLogOut << lorig.mNHCuts << " " << lorig.mNVCuts << " " << bestcost << endl;
				}
                l = lorig;
				if ((mcost = moveVCost(i, -d, maxactive)) < bestcost) {
					l.moveVCut(i, -d);
					best = l;
					bestcost = mcost;
					hit = true;
                    if (bestcost < 1E10)
                        mLogOut << lorig.mNHCuts << " " << lorig.mNVCuts << " " << bestcost << endl;
				}
			}

			// move each h cut by +/- d from orig, keep track of best cost
			// repeat for v
			// choose the one with lowest cost and return true
			// or false if none better

            if (hit)
                lorig = best;
            return hit;
        }

        void findBestLayoutNear(Layout &l, int dmax, int nactive) {
			int d;

			setPartitions(l.mNHCuts, l.mNVCuts, l.mHCuts, l.mVCuts);

            bool mademove;
			do {
				mademove = false;
				d = dmax;
				do {
					while(improvePartition(d, l, nactive)) {  // this should keep making individual h/v cut moves by +/- d until no improvement
						mademove = true;
					}
					d /= 2;
				} while (d>0);
			} while (mademove);
        }


	// can have many partitions, but never have nactive > target
	void findBestXY(int nx, int ny, int maxactive, Layout &l) {
		// int ntarget = maxactive;
		Layout best;


		int nh = nx;
		int nv = ny;

                findApproxLayout(nh, nv, l); 
                findBestLayoutNear(l, 32, maxactive);

                //cout << endl;
                //cout << "Final Best" << endl;
				//cout << l.mMaxCost << " " << l.mNHCuts << " " << l.mNVCuts << " " << l.mNActivePartitions << endl;
				//for (int i=0; i<nh; i++) cout << l.mHCuts[i]+MINTHICKX << " ";
				//cout << endl;
				//for (int i=0; i<nv; i++) cout << l.mVCuts[i]+MINTHICKY << " ";
				//cout << endl;
				//cout << endl;
                                
                // cout << l.mNHCuts << " " << l.mNVCuts << " " << l.mNActivePartitions << endl;
		//	}
		//}
		// cout << best.mNHCuts << " " << best.mNVCuts << " " << best.mNActivePartitions << endl;
	}


    void dump(char *fstem, Layout &lbest, int ntarget) {
        findCost(lbest);
        for (int i=0; i<lbest.mNPartitions; i++)
            lbest.mCosts[i] = mPartitions[i].mCost;

        int maxgwidth = 0;
        int maxgheight = 0;
        for (int i=0; i<lbest.mNHCuts; i++)
            if (maxgwidth<lbest.mHCuts[i])
                maxgwidth = lbest.mHCuts[i];

        for (int i=0; i<lbest.mNVCuts; i++)
            if (maxgheight<lbest.mVCuts[i])
                maxgheight = lbest.mVCuts[i];

        int gzwidth = 4;
        ofstream lorp;

        char buff[1024];

        sprintf(buff, "LORP_%s.INP", fstem);
        lorp.open(buff);

        lorp << buff << endl;
        lorp << " LORP input data used for decomposing domain into optimal configuration with " << ntarget << " processors" << endl;
        lorp << " NH NV NActive Cost TotalAreaWeight TileAreaWeight TileGhostWeight TileWaterWeight TileNearLandWeight" << endl;
        lorp << " " << lbest.mNHCuts << " " << lbest.mNVCuts << " " << lbest.mNActivePartitions << " " << lbest.mMaxCost << " " << mTotalAreaWeight << " " << 
                mTileAreaWeight << " " << mTileGhostWeight << " " << mTileWaterWeight << " " << mTileNearLandWeight << endl;
        cout << lbest.mNHCuts << " " << lbest.mNVCuts << " " << lbest.mNActivePartitions << " " << lbest.mMaxCost << " " << mTotalAreaWeight << " " << 
                mTileAreaWeight << " " << mTileGhostWeight << " " << mTileWaterWeight << " " << mTileNearLandWeight << endl;


        lorp << "*   IC:    an array of NPARTX size setting x-width for each partition" << endl;

        for (int i=0; i<lbest.mNHCuts-1; i++)
            lorp << lbest.mHCuts[i]+MINTHICKX+gzwidth << endl;
        lorp << lbest.mHCuts[lbest.mNHCuts-1]+gzwidth << endl;

        lorp << "*   JC:    an array of NPARTY size setting y-height for each partition" << endl; 
        for (int i=0; i<lbest.mNVCuts-1; i++)
            lorp << lbest.mVCuts[i]+MINTHICKY+gzwidth << endl;
        lorp << lbest.mVCuts[lbest.mNVCuts-1]+gzwidth << endl;

        lorp << "*   List of active partition IDs, 0-based" << endl;
        for (int i=0; i<lbest.mNPartitions; i++) {
            if (lbest.mCosts[i]>0) {
                lorp << i << endl;
            }
        }

        lorp << "*   Lookup table to go from your partition id to node id (so has full N entries and some are -1, 0-based)" << endl;
        int node = 0;
        for (int i=0; i<lbest.mNPartitions; i++) {
            int outval = -1;
            if (lbest.mCosts[i]>0) {
                outval = node;
                node++;
            }
            lorp << outval << endl;
        }

        lorp.close();

        ofstream cellview;
        char cviewname[1024];
        sprintf(cviewname, "cview_%s.txt", fstem);
        cellview.open(cviewname);

        int m = mHeight-1;
        for (int j=lbest.mNVCuts-1; j>=0; j--) {
            for (int jj=lbest.mVCuts[j]-1; jj>=0; jj--) {
                int n = 0;
                for (int i=0; i<lbest.mNHCuts; i++) {
                    for (int ii=0; ii<lbest.mHCuts[i]; ii++) {
                        cellview.put(char(mType[m][n]&0xf)+'0');
                        n++;
                    }
                    cellview.put(' ');
                }
                cellview << endl;
                m--;
            }
            cellview << endl;
        }
        cellview.close();
    }




    void doIt(char *cellfname, char *fstem, int ntarget) {
        Layout l;
        Layout lbest;
        lbest.mMaxCost = 1E90;
        char lstem[1024];

		loadTypes(cellfname);

        char buff[1024];

        sprintf(buff, "%s_%d_h.log", fstem, ntarget);
        sprintf(lstem, "%s_%d_h", fstem, ntarget);
        mLogOut.open(buff);

        ofstream clog;
        sprintf(buff, "%s_cost_%d_h.log", fstem, ntarget);
        clog.open(buff);

        for (int nv=1; nv<=1; nv++) {
            for (int nh=ntarget; nh<=ntarget; nh++) {
                if (nv*nh<ntarget)
                    continue;
                findBestXY(nh, nv, ntarget, l);
                if (l.mMaxCost < lbest.mMaxCost)
                    lbest = l;
            }
        }


        clog << endl;
        clog << " Tile Costs" << endl;
        for (int j=lbest.mNVCuts-1; j>=0; j--) {
            int npart = j*lbest.mNHCuts;
            clog << " ";
            for (int i=0; i<lbest.mNHCuts; i++) {
                clog << lbest.mCosts[npart] << " ";
                npart++;
            }
            clog << endl;
        }
        clog << endl;
        clog.close();

        mLogOut.close();

        dump(lstem, lbest, ntarget);

/////////////////////////////////////////////////

        lbest.mMaxCost = 1E90;


        sprintf(buff, "%s_%d_v.log", fstem, ntarget);
        mLogOut.open(buff);

        sprintf(buff, "%s_cost_%d_v.log", fstem, ntarget);
        sprintf(lstem, "%s_%d_v", fstem, ntarget);
        clog.open(buff);

        for (int nv=ntarget; nv<=ntarget; nv++) {
            for (int nh=1; nh<=1; nh++) {
                if (nv*nh<ntarget)
                    continue;
                findBestXY(nh, nv, ntarget, l);
                if (l.mMaxCost < lbest.mMaxCost)
                    lbest = l;
            }
        }

        clog << endl;
        clog << " Tile Costs" << endl;
        for (int j=lbest.mNVCuts-1; j>=0; j--) {
            int npart = j*lbest.mNHCuts;
            clog << " ";
            for (int i=0; i<lbest.mNHCuts; i++) {
                clog << lbest.mCosts[npart] << " ";
                npart++;
            }
            clog << endl;
        }
        clog << endl;
        clog.close();

        mLogOut.close();

        dump(lstem, lbest, ntarget);

/////////////////////////////////////////

        lbest.mMaxCost = 1E90;

        sprintf(buff, "%s_%d_r.log", fstem, ntarget);
        mLogOut.open(buff);

        sprintf(buff, "%s_cost_%d_r.log", fstem, ntarget);
        sprintf(lstem, "%s_%d_r", fstem, ntarget);
        clog.open(buff);

        for (int nv=1; nv<=(ntarget+2); nv++) {
            for (int nh=(ntarget/nv); nh<=(ntarget+2); nh++) {
                if (nv*nh<ntarget)
                    continue;
                findBestXY(nh, nv, ntarget, l);
                clog << nh << " " << nv << " " << l.mMaxCost << endl;
                if (l.mMaxCost < lbest.mMaxCost)
                    lbest = l;
            }
        }

        clog << endl;
        clog << " Tile Costs" << endl;
        for (int j=lbest.mNVCuts-1; j>=0; j--) {
            int npart = j*lbest.mNHCuts;
            clog << " ";
            for (int i=0; i<lbest.mNHCuts; i++) {
                clog << lbest.mCosts[npart] << " ";
                npart++;
            }
            clog << endl;
        }
        clog << endl;
        clog.close();

        mLogOut.close();

        dump(lstem, lbest, ntarget);
    }

};
