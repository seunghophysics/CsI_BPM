//#ifndef __CINT__

/* This script converts .pcoraw 2 by 2 pixel mode file to .root one.
 * The last contains a TTree object named "tree".  It has two branches.
 * t1 --- a TH2S * histogram,  a data-time stamp is situated in the left
 * top corner of the hisogram.  picnum --- an integer variable maches
 * the histogram counting number drawn on the pitcuture itself.
 *
 * Usage:
 * type in terminal:
 * $ root 'pco2root2x2.C+("/in/put/file.pcoraw", "/out/put/file.root")' -b -q
 * NOTE  Be shure that right file format are specified.
 *
 * The code was tested with ROOT 6.06/00.
 *
 * History of modification:
 * An initial code taken from a pco example.
 * Modified by BongHo K.
 * Modified by Gosha R.
 *
 * 2016-03-08
 */

#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <time.h>
#include <math.h>
// #include <chrono>
// #include <ctime>

#include "pco2root2x2.h"// contains definitions of structs.

using namespace std;
#define HEADERLEN 128
//#define BASEFILE
// #define PIXELNUM = 480000
const long PIXELNUM = 1920000;
const long PIXELNUM2 = 3840000;
const int HISTW = 1600;
const int HISTH = 1200;
const int PICNUMSTEP = 1194;
/*
*/
/*
const long PIXELNUM = 480000;
const long PIXELNUM2 = 960000;
const int HISTW = 800;
const int HISTH = 600;
const int PICNUMSTEP = 594;
*/

int lastpicnum = -1;


int getpicnum(TH2S * t1){
    int picnum = 0;
    for(int q=0; q<7; q++){
        if(t1->GetBinContent(67-8*q, PICNUMSTEP+3)>0 &&
                t1->GetBinContent(67+5-8*q, PICNUMSTEP+3)>0)
            picnum += pow(10, q) * 0;
        if(t1->GetBinContent(67+4-8*q, PICNUMSTEP+6)==0)
            picnum += pow(10, q) * 1;
        if(t1->GetBinContent(67+2-8*q, PICNUMSTEP+2)>0 &&
                t1->GetBinContent(67+3-8*q, PICNUMSTEP+2)>0 &&
                t1->GetBinContent(67+4-8*q, PICNUMSTEP+2)==0)
            picnum += pow(10, q) * 2;
        if(t1->GetBinContent(67+1-8*q, PICNUMSTEP+3)==0 &&
                t1->GetBinContent(67+2-8*q, PICNUMSTEP+3)>0)
            picnum += pow(10, q) * 3;
        if(t1->GetBinContent(67+3-8*q, PICNUMSTEP)==0 &&
                t1->GetBinContent(67+4-8*q, PICNUMSTEP)>0)
            picnum += pow(10, q) * 4;
        if(t1->GetBinContent(67+5-8*q, PICNUMSTEP+6)>0 &&
                t1->GetBinContent(67+5-8*q, PICNUMSTEP+5)==0)
            picnum += pow(10, q) * 5;
        if(t1->GetBinContent(67+5-8*q, PICNUMSTEP+4)==0 &&
                t1->GetBinContent(67+5-8*q, PICNUMSTEP+5)>1 &&
                t1->GetBinContent(67+5-8*q, PICNUMSTEP+6)==0)
            picnum += pow(10, q) * 6;
        if(t1->GetBinContent(67+1-8*q, PICNUMSTEP)>1 &&
                t1->GetBinContent(67+2-8*q, PICNUMSTEP)==0)
            picnum += pow(10, q) * 7;
        if(t1->GetBinContent(67-8*q, PICNUMSTEP+2)>0 &&
                t1->GetBinContent(67-8*q, PICNUMSTEP+3)==0 &&
                t1->GetBinContent(67-8*q, PICNUMSTEP+4)>0)
            picnum += pow(10, q) * 8;
        if(t1->GetBinContent(67-8*q, PICNUMSTEP+2)==0 &&
                t1->GetBinContent(67-8*q, PICNUMSTEP+3)==0 &&
                t1->GetBinContent(67-8*q, PICNUMSTEP+4)>0)
            picnum += pow(10, q) * 9;
    }
    cout << "The picture number is " << picnum << ".\n";
    return picnum;
}


// int pictureread(ifstream *pcofile,
void pictureread(ifstream *pcofile,
        unsigned short *picdata,
        int nulimit,
        int shiftv, 
        TTree *rtree,
        TH2S* t1,
        int * treepicnum,
        int save,
        int filetype,
        TH2S* baseh1,
        int basev){
    // std::cout << "nulimit = " << nulimit << ": shiftv = " << shiftv << endl;
    int binn;
    int picnum = -1;
    unsigned short picdatav;
    double testerr;
    int rtrn_fill;
    for(int nu=0; nu<nulimit; nu++){
        pcofile->read((char*)picdata, PIXELNUM2);
        if(save ==1 ){
            t1->Reset();
            for(int i=1; i<=HISTH; i++){
                for(int j=1; j<=HISTW; j++){
                    binn = (i-1)*HISTW + (j-1);
                    if(filetype == 0){
                        picdatav = picdata[binn]/4;
                        //cout<<picdatav<<endl;
                        // t1->SetBinContent(j, i, picdatav);
                        t1->SetBinContent(j, HISTH+1-i, picdatav);
                        if(basev == 1){
                            testerr = baseh1->GetBinError(j, i);
                            // t1->SetBinError(j, i, testerr);
                            t1->SetBinError(j, HISTH+1-i, testerr);
                        }
                    }
                    else
                        // t1->SetBinContent(j, i, picdata[binn]);
                        t1->SetBinContent(j, HISTH+1-i, picdata[binn]);
                }
            }
            picnum = getpicnum(t1);
            if(picnum == lastpicnum && false){
                break;// skip this and go to the *else* statement.
            }else{
                lastpicnum = picnum;
                (* treepicnum) = picnum;
                rtrn_fill = rtree->Fill();
                cout << rtrn_fill << " bytes were pushed in the tree.\n";
                cout << rtree->GetEntries() << " entries in the tree.\n";
            }
        }
        //		cout<<"function test = "<<pcofile->tellg()<<endl;          
        pcofile->seekg(shiftv, ios::cur);
    }
}


//int pco2root2x2(string input, string output){
 int main(int argc, char* argv[]){
    string input = argv[1];
    string output = argv[2];
    cout << "from = " << input.c_str() << endl;
    cout << "making = " << output.c_str() << endl;
    int save = 1; /////picture save yes :1 , no :0//////
    int filetype = 0;////will be decided at IFDV  filetype :1 (new), :0 (old)//////
    cout << "start = " << clock()/1000 << endl;
    long stime = time(NULL);
    clock_t clockstart = clock();// To measure processor time.
    time_t timestart = time(0);// To measure real time.

    ////////////////////20150619//////////////////
    TH2S *baseh1 = new TH2S("baseh1", "", HISTW, 0, HISTW, HISTH, 0, HISTH);
    int basev = 0;
#ifdef BASEFILE
        basev = 1;
        TFile *basefile = new TFile("../baseline_check/no_hv_withMylar_100ms_40MHz_2ADC_000_20150618_baseline.root");
        TTree *baset = (TTree*) basefile->Get("tree");
        baset->SetBranchAddress("normalt1", &baseh1);
        baset->GetEntry(0);
#endif
    /////////////////////////////////////////////

    TFile *rfile = new TFile(output.c_str(), "RECREATE");
    TTree *rtree = new TTree("tree", "tree");
    TH2S *t1 = new TH2S("t1", "t1", HISTW, 0, HISTW, HISTH, 0, HISTH);
    int treepicnum;
    // rtree->Branch("t1", &t1, 32000, 0);
    rtree->Branch("t1", &t1);
    rtree->Branch("picnum", &treepicnum, "picnum/I");

    ifstream pcofile;
    pcofile.open(input.c_str(), ios::in | ios::binary);

    pcofile.seekg(0, ios::end);
    long long pcofilesize = (long long)pcofile.tellg();
    // cout << "PCO file size is " << pcofilesize << "\n";
    // cout << "The size of non-pictures is " << (pcofilesize-(2*PIXELNUM+908)*10000) / 10000. << "\n";
    long picturenum = (long)(pcofile.tellg())/(2*PIXELNUM+908);
    cout << "total picture number is = " << picturenum << endl;
    pcofile.seekg(0, ios::beg);
    // cout << "//////////////////first header///////////////" << endl;
    unsigned short ttype;// = TIFFEntry.ftype;
    pcofile.read((char*)&ttype, sizeof(ttype));
    // cout << "TIFFEntry.ftype " << ttype << endl;
    unsigned short ttfield;// = TIFFEntry.TagField;
    pcofile.read((char*)&ttfield, sizeof(ttfield));
    // cout << "TIFFEntry.TagField " << ttfield << endl;
    unsigned short tlen;// = TIFFEntry.length;
    pcofile.read((char*)&tlen, sizeof(tlen));
    // cout << "TIFFEntry.length " << tlen << endl;//":"<<lh<<endl;
    unsigned short toff;// = TIFFEntry.offset;
    pcofile.read((char*)&toff, sizeof(toff));
    // cout << "TIFFEntry.offset " << toff << endl;
    unsigned long ifdv;
    pcofile.read((char*)&ifdv, 8);//sizeof(ifdv));
    // cout << "tag offset point = " << ifdv << endl;
    // cout << "file type = " << filetype << endl;
    long firstoff = pcofile.tellg();

    pcofile.seekg(ifdv, ios::beg);
    unsigned long taglength;
    pcofile.read((char*)&taglength, sizeof(long));
    // cout << "length of first tag : " << taglength << endl;
    // long secondimage = ifdv + (taglength*20) + (8*2);// Unused.
    unsigned short tfv, ftv;
    unsigned long ltv, osv;
    for(int i=0; i<taglength; i++){
        pcofile.read((char*)&tfv, 2);
        pcofile.read((char*)&ftv, 2);
        pcofile.read((char*)&ltv, 8);
        pcofile.read((char*)&osv, 8);
    }
    // cout << "///////////////////////////////////" << endl;
    pcofile.seekg(firstoff, ios::beg);
    unsigned char *testt = new unsigned char[49];
    pcofile.read((char*)testt, 49);//sizeof(testt));
    // int count = 0;// Unused.
    // int countn = 0;// Unused.
    /*for(int i=0; i<49; i++){
        cout << testt[i];//<<endl;
    }
    cout << endl;*/
    // cout << "data position0 = " << pcofile.tellg() << endl;

    Bild testb;

    unsigned short *tests = new unsigned short[PIXELNUM];
    pictureread(&pcofile, tests, 1, 0, rtree, t1, &treepicnum, save, filetype, baseh1, basev);

    // picturenum = 1;// Just fot the test.
    if(picturenum == 1){
        cout << "The number of pictures is 1 and we finish here.\n";
        rtree->Write();
        rfile->Write();
        rfile->Close();
        cout << "end = " << clock()/1000 << endl;
        cout << time(NULL)-stime << endl;
        return 0;
    }
    pcofile.seekg(4, ios::cur);
    pcofile.read((char*)&testb, sizeof(Bild));
    int shiftbv = 684-(2*PIXELNUM+969)+ifdv;//1164 - 3840973 + ifdv;
    pcofile.seekg(shiftbv, ios::cur);
    ///////////////second tag///////////////////////
    pcofile.read((char*)&taglength, 8);
    // cout << "The tag length is " << taglength << ".\n";
    // unsigned short tfv, ftv;
    // unsigned long ltv, osv;
    for(int i=0; i<taglength; i++){
        pcofile.seekg(20, ios::cur);
    }
    // cout << "///////////////////////////////////" << endl;
    // cout << "here is = " << pcofile.tellg() << endl;
    pcofile.seekg(8, ios::cur);
    // cout << "data position1 = " << pcofile.tellg() << endl;
    unsigned short *tests1 = new unsigned short[PIXELNUM];
    pictureread(&pcofile, tests1, 1, 0, rtree, t1, &treepicnum, save, filetype, baseh1, basev);
    cout << rtree->GetEntries() << " entries in the tree.\n";
    if(picturenum == 2){
        cout << "The number of pictures is 2 and we finish here.\n";
        rtree->Write();
        rfile->Write();
        rfile->Close();
        cout << "end = " << clock()/1000 << endl;
        cout << time(NULL)-stime << endl;
        return 0;
    }

    ////////////////caseb/////////////////////////////////////
    int shiftv = 904-(2*PIXELNUM+969)+ifdv;
    pcofile.seekg(shiftv, ios::cur);
    ////////////////////////////////////////////////////
    // cout << "data position2 = " << pcofile.tellg() << endl;


    pictureread(&pcofile, tests1, 1, 0, rtree, t1, &treepicnum, save, filetype, baseh1, basev);
    cout << rtree->GetEntries() << " entries in the tree.\n";
    if(picturenum == 3){
        cout << "The number of pictures is 3 and we finish here.\n";
        rtree->Write();
        rfile->Write();
        rfile->Close();
        cout << "end = " << clock()/1000 << endl;
        cout << time(NULL)-stime << endl;
        return 0;
    }

    //break;
    //////////picture reading from pic2 to pic99/////////////////
    unsigned short *tests2 = new unsigned short[PIXELNUM];
    long picturenumt = picturenum/100;
    int nulimit;
    if (picturenumt >=1)
        nulimit = 97;//98; 
    else 
        nulimit = picturenum-3;
    for(int nu = 0; nu<nulimit; nu++){
        pcofile.seekg(shiftv, ios::cur);
        // cout << "data position " << nu+3 << " = " << pcofile.tellg() << endl;
        pictureread(&pcofile, tests2, 1, 0, rtree, t1, &treepicnum, save, filetype, baseh1, basev);
        cout << rtree->GetEntries() << " entries in the tree.\n";
    }

    // picturenumt = 0;// Added specialy to test with 2x2
    if(picturenumt == 0){
        cout << "The number of hundreds of pictures is 0 and we finish here.\n";
        cout << "Saving :  " << rtree->GetEntries() << " entries in the tree.\n";
        int rtrn_tree_write = rtree->Write();
        cout << "Saving :  " << rtree->GetEntries() << " entries in the tree.\n";
        cout << rtrn_tree_write << " bytes were written in the directory/file.\n";
        rfile->Write();
        rfile->Close();
        clock_t clockstop = clock();// To measure processor time.
        time_t timestop = time(0);// To measure real time.
        double timecpu = double(clockstop-clockstart) / CLOCKS_PER_SEC;
        std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
                  // << 1000.0 * (clockstop-clockstart) / CLOCKS_PER_SEC << " ms "
                  << timecpu << " s "
                  << "( " << 100.0/timecpu << " Hz )\n"
                  << "Wall clock time passed: "
                  << (timestop-timestart)*1000.0
                  << " ms "
                  << "( " << 100.0 / (timestop-timestart) << " Hz )\n";
        //// Try to solve the problem with tree saving only one entry. ////
        /*
        TFile Fout("data/20170808/converted/tt.root", "recreate");
        rtree->Write("tt");
        Fout.Close();
        */
        return 0;
    }


    //break;
    ////////////////casea///////////////////////////////////////
    // cout << "before reading bild = " << pcofile.tellg() << endl;
    //break;
    pcofile.seekg(4, ios::cur);

    pcofile.read((char*)&testb, sizeof(Bild));
    // cout << "year = " << testb.sTime.wYear << endl;
    // cout << testb.iXRes << ":" << testb.iYRes << endl;
    // cout << "align = 1 = " << testb.bAlignUpper << endl;
    unsigned long gb = testb.dGammaLut;
    // cout << gb << ":" << testb.dGammaLutC << ":" << testb.dGammaLut2 << ":" << testb.dGammaLutC2 << endl;
    // cout << "bitres = " << testb.iBitRes << endl;
    ////break;

    //////////tag info///////////////
    pcofile.seekg(shiftbv, ios::cur);
    // cout << "tag loading ....................from " << pcofile.tellg() << endl;
    int addshift = -1980;
    pcofile.seekg(addshift, ios::cur);
    int totaltag;
    totaltag = taglength*20;
    for(int tagn=0; tagn<99; tagn++){
        pcofile.read((char*)&taglength, 8);
        pcofile.seekg(totaltag, ios::cur);
        pcofile.seekg(8, ios::cur);
        //cout<<pcofile.tellg()<<endl;
    }
    // cout << "taglength = " << taglength << "\n";
    // cout << "data position" << nulimit+3 << " = " << pcofile.tellg() << endl;
    //break;
    //////////////////////////////
    int curr = 0;
    int allrun = 0;
    // picturenum -= 2;
    int leftpict = picturenum % 100;
    nulimit = 100;
    // cout << "number of pic = " << nulimit << " : " << leftpict << endl;
    // int addshift = -1980;
    // pcofile.seekg(addshift, ios::cur);
    // cout << "The additional shift is addshift = " << addshift << "\n";
    int escaper = 4;
    do
    {
        /*
        if(escaper == 0)
            break;
        escaper -= 1;
        cout << "escaper = " << escaper  << "\n";
        */

        if(curr == (picturenumt-1)){
            // cout << "almost done" << endl;
            nulimit = leftpict;
            //allrun = 1;
        }

        // cout << "Before pictureread the position is " << pcofile.tellg() << ".\n";
        pictureread(&pcofile, tests2, nulimit, shiftv, rtree, t1, &treepicnum, save, filetype, baseh1, basev);
        // pcofile.seekg(shiftv,ios::cur);
        // cout << "After pictureread the position is " << pcofile.tellg() << ".\n";
        pcofile.seekg(addshift-20, ios::cur);

        if(curr == (picturenumt-1) && nulimit == leftpict)
            allrun = 1;

        int totaltag = taglength*20;
        totaltag = 12*20;
        // cout << "totaltag = " << totaltag << "\n";
        for(int tagn = 0; tagn<nulimit; tagn++){
            pcofile.read((char*)&taglength, 8);
            pcofile.seekg(totaltag, ios::cur);
            pcofile.seekg(8, ios::cur);
        }
        curr++;
        cout << "curr = " << curr << "\n";
        // cout << "After tag reading the position is " << pcofile.tellg() << ".\n";
    }while (!(allrun = 1 && nulimit == leftpict));

    rtree->Write();
    rfile->Write();
    rfile->Close();
    cout << "end = " << clock()/1000. << endl;
    cout << time(NULL)-stime << endl;
    clock_t clockstop = clock();// To measure processor time.
    time_t timestop = time(0);// To measure real time.
    double timecpu = double(clockstop-clockstart) / CLOCKS_PER_SEC;
    std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
              // << 1000.0 * (clockstop-clockstart) / CLOCKS_PER_SEC << " ms "
              << timecpu << " s "
              << "( " << 100.0/timecpu << " Hz )\n"
              << "Wall clock time passed: "
              << (timestop-timestart)*1000.0
              << " ms "
              << "( " << 100.0 / (timestop-timestart) << " Hz )\n";
    return 0;
}

