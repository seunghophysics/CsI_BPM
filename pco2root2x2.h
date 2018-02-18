#define HEADERLEN 128

typedef struct
{
  char ucPco[4];
  unsigned int uiFileLen;
  unsigned int uiHeaderLen;
  unsigned int uiXRes;
  unsigned int uiYRes;
  unsigned int uiLutSign;
  unsigned int uiColor;
  unsigned int uiBWMin;
  unsigned int uiBWMax;
  unsigned int uiBWLut;
  unsigned int uiRMin;
  unsigned int uiRMax;
  unsigned int uiGMin;
  unsigned int uiGMax;
  unsigned int uiBMin;
  unsigned int uiBMax;
  unsigned int uiColLut;
  unsigned int uiDS;
  unsigned int uiDummy[HEADERLEN-18];
}B16_HEADER;
typedef unsigned long DWORD;
typedef unsigned short WORD;

typedef struct _SYSTEMTIME {
    WORD wYear;
    WORD wMonth;
    WORD wDayOfWeek;
    WORD wDay;
    WORD wHour;
    WORD wMinute;
    WORD wSecond;
    WORD wMilliseconds;
} SYSTEMTIME, *PSYSTEMTIME, *LPSYSTEMTIME;

typedef struct
{
//    WORD *pic12 = new unsigned short[testsize];
//  WORD *pic12;                         // original image pointer
  SYSTEMTIME sTime;                    // shows the exact time stamp of the image  // 20 byte
  int        iTicks;                   // milliseconds gone after start in stime
  int        iXRes;                    // X Resolution
  int        iYRes;                    // Y Resolution                             // 32 byte
  char cText[40];                      // text which should be placed inside the image// 72 byte
  bool       bDouble;                  // shows if pic is Doubleshutter image
  bool       bDummy[3];                // since bool is only one byte, we've to fill up to four bytes// 76 byte
  int        iBWMin;                   // Lut bw min                               // 80 byte
  int        iBWMax;                   // Lut bw max
  int        iBWLut;                   // Lut lin log
  int        iRMin;                    // red min                                  // 92 byte
  int        iRMax;                    // red max
  int        iGMin;                    // green min                                // 100 byte
  int        iGMax;                    // green max
  int        iBMin;                    // blue min
  int        iBMax;                    // blue max                                 // 112 byte
  int        iColLut;                  // Lut lin log color
  int        iColor;                   // image from Color-CCD: 1 otherwise 0      // 120 byte
  int        iVersion;                 // Version of b16 extended info
  int        iBWMin2;                   // Lut bw min
  int        iBWMax2;                   // Lut bw max                              // 132 byte
  int        iBWLut2;                   // Lut lin log
  int        iRMin2;                    // red min                                 // 140 byte
  int        iRMax2;                    // red max
  int        iGMin2;                    // green min
  int        iGMax2;                    // green max                               // 152 byte
  int        iBMin2;                    // blue min
  int        iBMax2;                    // blue max                                // 160 byte
  int        iColLut2;                  // Lut lin log color
  bool       bAlignUpper;               // Align MSB (0-MSB is bit14, 1-MSB is bit 16)
  bool       bDummy2[3];                // since bool is only one byte, we've to fill up to four bytes // 168 byte
//  double
  double     dGammaLut;                 // Gamma value b/w
  double     dGammaLutC;                // Gamma value color
  double     dGammaLut2;                // Gamma value b/w 2
  double     dGammaLutC2;               // Gamma value color 2                     // 200 byte
  int        iColorPatternType;         // Demosaicking type for the color pattern
  int        iBitRes;                   // Bit resolution of image                 // 208 byte
  double     dSaturation;               // Color saturation common for both ds images // 216 byte
}Bild;// ACHTUNG: noch 172 Bytes frei, sonst muss headerlen in file12 angepasst werden!

typedef struct
{
    unsigned short TagField;
    unsigned short ftype;
    unsigned short length;
    unsigned short offset;
}TE;

