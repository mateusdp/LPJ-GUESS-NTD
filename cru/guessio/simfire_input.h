//////////////////////////////////////////////////////////////////////////////////////
// SIMFIREINPUT.H
// Header file for input from a fast data archive
// Created automatically by FastArchive on Tue Jun 16 13:43:53 2015
//
// For a description of the data, please see getsimfiredata (modules/simfire.cpp)
//
// The following #includes should appear in your source code file:
//
//   #include <stdio.h>
//   #include <stdlib.h>
//   #include <string.h>
//   #include "SimfireInput.h"
//
// Functionality to retrieve data from the archive is provided by class SimfireInputArchive.
// The following public functions are provided:
//
// bool open(char* filename)
//   Attempts to open the specified file as a fast data archive. The format must be
//   exactly compatible with this version of SimfireInput.h (normally the archive and
//   header file should have been produced together by the same program using class
//   CFastArchive). Returns false if the file could not be opened or had format
//   errors. open() with no argument is equivalent to open("SimfireInput.bin").
//
// void close()
//   Closes the archive (if open).
//
// bool rewind()
//   Sets the file pointer to the first record in the archive file. Returns false if
//   no archive file is currently open.
//
// bool getnext(SimfireInput& obj)
//   Retrieves the next record in the archive file and advances the file pointer to
//   the next record. Data are written to the member variables of obj. Returns false if
//   no archive file is currently open or if the file pointer is beyond the last
//   record. Use rewind() and getnext() to retrieve data sequentially from the archive.
//
// bool getindex(SimfireInput& obj)
//   Searches the archive for a record matching the values specified for the index
//   items (lon and lat) in obj. If a matching record is found, the data are
//   written to the member variables of obj. Returns true if the archive was open and
//   a matching record was found, otherwise false. The search is iterative and fast.
//
// Sample program:
//
//   SimfireInputArchive ark;
//   SimfireInput data;
//   bool success,flag;
//
//   // Retrieve all records in sequence and print values of lon and lat:
//
//   success=ark.open("SimfireInput.bin");
//   if (success) {
//      flag=ark.rewind();
//      while (flag) {
//         flag=ark.getnext(data);
//         if (flag)
//            printf("Loaded record: lon=%g, lat=%g\n",data.lon,data.lat);
//      }
//   }
//   
//   // Look for a record with lon=-179.75, lat=-89.75:
//
//   data.lon=-179.75;
//   data.lat=-89.75;
//   success=ark.getindex(data);
//   if (success) printf("Found it!\n");
//   else printf("Not found\n");
//
//   ark.close();


struct SimfireInput {

	// Index part

	double lon;
	double lat;

	// Data part

	double igbp_class[1];
	double pop_density[57];
	double monthly_burned_area[12];
};


const long SIMFIREINPUT_NRECORD=259200;
const int SIMFIREINPUT_DATA_LENGTH=258;
const int SIMFIREINPUT_INDEX_LENGTH=8;
const int SIMFIREINPUT_HEADERSIZE=563;
unsigned char SIMFIREINPUT_HEADER[SIMFIREINPUT_HEADERSIZE-4]={
	0x01,0x02,0x33,0x00,0x00,0x00,0x08,0x00,0x02,0x04,0x6C,0x6F,0x6E,0x00,0x2D,0x31,0x37,0x39,0x2E,0x37,
	0x35,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,
	0x00,0x00,0x00,0x00,0x00,0x00,0x31,0x37,0x39,0x2E,0x37,0x35,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,
	0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x30,0x2E,
	0x30,0x35,0x00,0x35,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,
	0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x0D,0x04,0x6C,0x61,0x74,0x00,0x2D,0x38,0x39,0x2E,
	0x37,0x35,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,
	0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x38,0x39,0x2E,0x37,0x35,0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,
	0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
	0x30,0x2E,0x30,0x35,0x00,0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,
	0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x0C,0x00,0x00,0x01,0x02,0x00,0x03,0x0B,
	0x69,0x67,0x62,0x70,0x5F,0x63,0x6C,0x61,0x73,0x73,0x00,0x30,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x5D,
	0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,
	0x00,0x00,0x00,0x31,0x37,0x00,0x00,0x00,0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,
	0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x31,0x00,0x00,0x00,0x00,
	0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,
	0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x05,0x00,0x00,0x00,0x01,0x0C,0x70,0x6F,0x70,0x5F,0x64,0x65,0x6E,
	0x73,0x69,0x74,0x79,0x00,0x30,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,
	0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x32,0x30,0x30,
	0x39,0x37,0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,
	0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x31,0x65,0x2D,0x30,0x35,0x00,0x00,0x00,0x5D,0xBF,0x53,
	0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,
	0x00,0x1F,0x00,0x00,0x00,0x39,0x0B,0x6D,0x6F,0x6E,0x74,0x68,0x6C,0x79,0x5F,0x62,0x61,0x00,0x30,0x00,
	0x00,0x00,0x00,0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,
	0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x31,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x5D,0xBF,
	0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,
	0x00,0x00,0x31,0x65,0x2D,0x30,0x37,0x00,0x00,0x00,0x5D,0xBF,0x53,0xEE,0xE8,0x7F,0x00,0x00,0x10,0x40,
	0xEB,0xAF,0xFF,0x7F,0x00,0x00,0xEC,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x18,0x00,0x00,0x00,0x0C};


class SimfireInputArchive {

private:

	FILE* pfile;
	long recno;
	long datano;
	unsigned char pindex[SIMFIREINPUT_INDEX_LENGTH];
	unsigned char pdata[SIMFIREINPUT_DATA_LENGTH];
	bool iseof;

	long readbin(int nbyte) {

		unsigned char buf[4];
		long mult[4]={0x1000000,0x10000,0x100,1},val;
		int i;

		fread(buf,nbyte,1,pfile);

		val=0;
		for (i=0;i<nbyte;i++) {
			val+=buf[i]*mult[4-nbyte+i];
		}

		return val;
	}

	void getindex(long n) {

		fseek(pfile,SIMFIREINPUT_INDEX_LENGTH*(n-recno),SEEK_CUR);
		fread(pindex,SIMFIREINPUT_INDEX_LENGTH,1,pfile);
		datano=pindex[SIMFIREINPUT_INDEX_LENGTH-4]*0x1000000+pindex[SIMFIREINPUT_INDEX_LENGTH-3]*0x10000+
			pindex[SIMFIREINPUT_INDEX_LENGTH-2]*0x100+pindex[SIMFIREINPUT_INDEX_LENGTH-1];
		recno=n+1;
		iseof=(recno==SIMFIREINPUT_NRECORD);
	}

	void getdata() {

		fseek(pfile,SIMFIREINPUT_INDEX_LENGTH*-recno+(datano-SIMFIREINPUT_NRECORD)*SIMFIREINPUT_DATA_LENGTH,SEEK_CUR);
		fread(pdata,SIMFIREINPUT_DATA_LENGTH,1,pfile);
		fseek(pfile,SIMFIREINPUT_DATA_LENGTH*(SIMFIREINPUT_NRECORD-datano-1)+SIMFIREINPUT_INDEX_LENGTH*recno,SEEK_CUR);
	}

	double popreal(unsigned char* bits,int nbyte,int nbit,double scalar,double offset) {

		unsigned char buf;
		int nb=nbit/8,i;
		double rval=0.0;
		long mult[4]={1,0x100,0x10000,0x1000000};

		for (i=0;i<4;i++) {
			if (i<nb) rval+=bits[nbyte-i-1]*mult[i];
			else if (i==nb) {
				buf=bits[nbyte-i-1]<<(8-nbit%8);
				buf>>=8-nbit%8;
				rval+=buf*mult[i];
			}
		}

		for (i=nbyte-1;i>=0;i--) {
			if (i>=nb)
				bits[i]=bits[i-nb];
			else
				bits[i]=0;
		}

		nb=nbit%8;

		for (i=nbyte-1;i>=0;i--) {
			bits[i]>>=nb;
			if (i>0) {
				buf=bits[i-1];
				buf<<=8-nb;
				bits[i]|=buf;
			}
		}

		rval=rval*scalar+offset;

		return rval;
	}

	void bitify(unsigned char buf[4],double fval,double offset,double scalar) {

		long ival = (long)((fval-offset)/scalar + 0.5);
		buf[0]=(unsigned char)(ival/0x1000000);
		ival-=buf[0]*0x1000000;
		buf[1]=(unsigned char)(ival/0x10000);
		ival-=buf[1]*0x10000;
		buf[2]=(unsigned char)(ival/0x100);
		ival-=buf[2]*0x100;
		buf[3]=(unsigned char)(ival);
	}

	void merge(unsigned char ptarget[4],unsigned char buf[4],int bits) {

		int nb=bits/8;
		int i;
		unsigned char nib;
		for (i=0;i<4;i++) {

			if (i<4-nb)
				ptarget[i]=ptarget[i+nb];
			else
				ptarget[i]=0;
		}
		nb=bits%8;
		for (i=0;i<4;i++) {
			ptarget[i]<<=nb;
			if (i<4-1) {
				nib=ptarget[i+1]>>(8-nb);
				ptarget[i]|=nib;
			}
		}

		nb=bits/8;
		if (bits%8) nb++;
		for (i=1;i<=nb;i++)
			ptarget[4-i]|=buf[4-i];
	}

	int compare_index(unsigned char* a,unsigned char* b) {

		int i;
		for (i=0;i<4;i++) {
			if (a[i]<b[i]) return -1;
			else if (a[i]>b[i]) return +1;
		}

		return 0;
	}

	bool initialise(const char* filename) {

		int i;
		unsigned char* pheader;

		if (pfile) fclose(pfile);
		pfile=fopen(filename,"rb");
		if (!pfile) {
			printf("Could not open %s for input\n",filename);
			return false;
		}

		pheader=new unsigned char[SIMFIREINPUT_HEADERSIZE-4];
		if (!pheader) {
			printf("Out of memory\n");
			fclose(pfile);
			pfile=NULL;
			return false;
		}
		::rewind(pfile);
		fread(pheader,SIMFIREINPUT_HEADERSIZE-4,1,pfile);
		for (i=0;i<SIMFIREINPUT_HEADERSIZE-4;i++) {
			if (pheader[i]!=SIMFIREINPUT_HEADER[i]) {
				printf("Format of %s incompatible with this version of simfire_input.h\n",filename);
				fclose(pfile);
				pfile=NULL;
				delete pheader;
				return false;
			}
		}
		delete[] pheader;

		::rewind(pfile);
		fseek(pfile,SIMFIREINPUT_HEADERSIZE+SIMFIREINPUT_DATA_LENGTH*SIMFIREINPUT_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

public:

	SimfireInputArchive() {
		pfile=NULL;
	}

	~SimfireInputArchive() {
		if (pfile) fclose(pfile);
	}

	bool open(const char* filename) {
		return initialise(filename);
	}

	bool open() {
		return open("simfire_input.bin");
	}

	void close() {
		if (pfile) {
			fclose(pfile);
			pfile=NULL;
		}
	}

	bool rewind() {

		if (!pfile) return false;

		::rewind(pfile);
		fseek(pfile,SIMFIREINPUT_HEADERSIZE+SIMFIREINPUT_DATA_LENGTH*SIMFIREINPUT_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

	bool getnext(SimfireInput& obj) {

		if (!pfile || iseof) return false;

		int i;

		getindex(recno);
		getdata();

		obj.lat=popreal(pindex,4,12,0.05,-89.75);
		obj.lon=popreal(pindex,4,13,0.05,-179.75);

		for (i=11;i>=0;i--) obj.monthly_burned_area[i]=popreal(pdata,SIMFIREINPUT_DATA_LENGTH,24,1e-07,0);
		for (i=56;i>=0;i--) obj.pop_density[i]=popreal(pdata,SIMFIREINPUT_DATA_LENGTH,31,1e-05,0);
		for (i=0;i>=0;i--) obj.igbp_class[i]=popreal(pdata,SIMFIREINPUT_DATA_LENGTH,5,1,0);

		return true;
	}

	bool getindex(SimfireInput& obj) {

		if (!SIMFIREINPUT_NRECORD || !pfile) return false;

		// else

		unsigned char ptarget[4]={0,0,0,0};
		unsigned char buf[4];
		bitify(buf,obj.lon,-179.75,0.05);
		merge(ptarget,buf,13);
		bitify(buf,obj.lat,-89.75,0.05);
		merge(ptarget,buf,12);

		long begin = 0;
		long end = SIMFIREINPUT_NRECORD;

		while (begin < end) {
			long middle = (begin+end)/2;

			getindex(middle);
			getdata();

			int c = compare_index(pindex, ptarget);

			if (c < 0) {
				begin = middle + 1;
			}
			else if (c > 0) {
				end = middle;
			}
			else {

				for (int i=11;i>=0;i--) obj.monthly_burned_area[i]=popreal(pdata,SIMFIREINPUT_DATA_LENGTH,24,1e-07,0);
				for (int i=56;i>=0;i--) obj.pop_density[i]=popreal(pdata,SIMFIREINPUT_DATA_LENGTH,31,1e-05,0);
				for (int i=0;i>=0;i--) obj.igbp_class[i]=popreal(pdata,SIMFIREINPUT_DATA_LENGTH,5,1,0);

				return true;
			}
		}

		return false;
	}
};
