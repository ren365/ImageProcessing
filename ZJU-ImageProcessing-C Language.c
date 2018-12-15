#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define BITMAPFILEHEADERLENGTH 14   // The bmp FileHeader length is 14
#define BM 19778                    // The ASCII code for BM


//================My file 1 start=============//
#define MAXSIZE 10000
unsigned char Y[MAXSIZE][MAXSIZE],output_Y[MAXSIZE][MAXSIZE];
unsigned char U[MAXSIZE][MAXSIZE],output_U[MAXSIZE][MAXSIZE];
unsigned char V[MAXSIZE][MAXSIZE],output_V[MAXSIZE][MAXSIZE];
void RGBtoYUV();
void YUVtoGray();
void YUVtoRGB();
void ChangeTheLuminance();
#define SourceFile "C:\\Users\\Administrator\\Desktop\\air.bmp"
#define OutputFile1 "C:\\Users\\Administrator\\Desktop\\air1.bmp"
#define OutputFile2 "C:\\Users\\Administrator\\Desktop\\air2.bmp"
#define OutputFile3 "C:\\Users\\Administrator\\Desktop\\air3.bmp"
#define OutputFile4 "C:\\Users\\Administrator\\Desktop\\air4.bmp"
#define OutputFile5 "C:\\Users\\Administrator\\Desktop\\air5.bmp"
#define OutputFile6 "C:\\Users\\Administrator\\Desktop\\air6.bmp"
#define OutputFile7 "C:\\Users\\Administrator\\Desktop\\air7.bmp"
//================My file 1 end=============//

//================My file 2 start=============//
void YUVtoBinarization();
#define elementline 3
void Dilation_Lining(char dowhat);
void Erosion();
void Opening_Closing(char dowhat);
void BetweenTwoOperationProcessing();
//================My file 2 end=============//

//================My file 3 end=============//
void VisibilityEnhancement();
void HistogramEqualization_r();
void HistogramEqualization_g();
void HistogramEqualization_b();
void HistogramEqualization_Y();
//================My file 3 end=============//

//================My file 4 start==========//
void Scaling(float ratioX,float ratioY);
void Rotation(double angle);
double Gaussian(int radial,double lamda);
void RBF_Gaussian(unsigned char temple[],int x,int y,int paramenter);
void Shearing(float dx,float dy);
void mirror(char character);
void translation(int sizeX,int sizeY,int sizeX1,int sizeY1,int x,int y);
void translation1(int sizeX,int sizeY,int sizeX1,int sizeY1,int x,int y);
//================My file 4 end=============//

//================My file 5 start==========//
void meanFiltering(int paramenter);
void MeanFilter(unsigned char temple[],int x,int y,int paramenter);
unsigned char Laplacian_filter_Y(int x,int y);
void LaplacianFilter_Y(double paramenter);
void Meanfilter_Y(int paramenter);
//================My file 5 end=============//

//================My file 6 start==========//
void BilateralFilter(int paramenter,double lamda_range,double lamda_color);
void RBF_Bilateral(unsigned char temple[],int x,int y,int paramenter,double lamda_range,double lamda_color);
//================My file 6 end=============//

//================Someone's work===========//
/* Test the file is bmp file or not */
void bmpFileTest(FILE* fpbmp);
/* To get the OffSet of header to data part */
void bmpHeaderPartLength(FILE* fpbmp);
/* To get the width and height of the bmp file */
void BmpWidthHeight(FILE* fpbmp);
//get r,g,b data
void bmpDataPart(FILE* fpbmp);
// output data to corresponding txt file
void bmpoutput(FILE *fpout);
//===========Someone's work end===========//

//===========Rewrite File ================//
void Initialization();
//open a picture
FILE *openfile(char filename[]);
//creat a picture
FILE *writefile(char filename[]);
//add head to file
void addHeadertofile(FILE *input,FILE *output);
 //===========Rewrite File End================//
 
 //============TEST FUNCTION=================//
void testInputYUV();
void testOutputYUV();
void testRGB();
void NoneProcess();
 //=============TEST FUNCTION END============//

unsigned int OffSet = 0;    // OffSet from Header part to Data Part
long width;          // The Width of the Data Part
long height;         // The Height of the Data Part
unsigned char r[MAXSIZE][MAXSIZE], output_r[MAXSIZE][MAXSIZE];
unsigned char g[MAXSIZE][MAXSIZE], output_g[MAXSIZE][MAXSIZE];
unsigned char b[MAXSIZE][MAXSIZE], output_b[MAXSIZE][MAXSIZE];
unsigned char records[MAXSIZE][MAXSIZE];

int main(int argc, char* argv[])
{
	/* Open bmp file */
	unsigned char *fp_temp;
	
	FILE *fpbmp;
	FILE *fpout;
	Initialization();
	
//===============first picture -- Bilateral Filter==========
	
	/*Standerd operations for file i/o*/
	fpbmp=openfile(SourceFile);
	bmpDataPart(fpbmp);
	
	fpout=writefile(OutputFile1);
	addHeadertofile(fpbmp,fpout);
	
	/*your operations for your picture*/
	
	bmpoutput(fpout);
	
	/*Standerd operations for close*/
	fclose(fpbmp);
	fclose(fpout);
		
	return 0;
}


void bmpoutput(FILE* fpout)
{
	long i, j = 0;
	long stride;
	unsigned char* pixout = NULL;

	stride = (24 * width + 31) / 8;
	stride = stride / 4 * 4;
	pixout = malloc(stride);

	fseek(fpout, OffSet, SEEK_SET);


	for (j = 0; j<height; j++)
	{
		for (i = 0; i<width; i++)
		{
			pixout[i * 3 + 2] = output_r[height - 1 - j][i];
			pixout[i * 3 + 1] = output_g[height - 1 - j][i];
			pixout[i * 3] = output_b[height - 1 - j][i];
		}
		fwrite(pixout, 1, stride, fpout);

	}
}


void bmpDataPart(FILE* fpbmp)
{
	int i, j = 0;
	int stride;
	unsigned char* pix = NULL;

	fseek(fpbmp, OffSet, SEEK_SET);
	stride = (24 * width + 31) / 8;
	stride = stride / 4 * 4;
	pix = malloc(stride);

	for (j = 0; j<height; j++)
	{
		fread(pix, 1, stride, fpbmp);

		for (i = 0; i<width; i++)
		{
			r[height - 1 - j][i] = pix[i * 3 + 2];
			g[height - 1 - j][i] = pix[i * 3 + 1];
			b[height - 1 - j][i] = pix[i * 3];


			output_r[height - 1 - j][i] = 255;
			output_g[height - 1 - j][i] = 255;
			output_b[height - 1 - j][i] = 255;
		}

	}
}

void bmpFileTest(FILE* fpbmp)
{
	unsigned short bfType = 0;

	fseek(fpbmp, 0L, SEEK_SET);//seek_set 起始位置
	fread(&bfType, sizeof(char), 2, fpbmp);
	if (BM != bfType)
	{
		printf("This file is not bmp file.!!!\n");
		exit(1);
	}
}


/* To get the OffSet of header to data part */
void bmpHeaderPartLength(FILE* fpbmp)
{
	fseek(fpbmp, 10L, SEEK_SET);
	fread(&OffSet, sizeof(char), 4, fpbmp);
	printf("The Header Part is of length %d.\n", OffSet);
}


/* To get the width and height of the bmp file */
void BmpWidthHeight(FILE* fpbmp)
{
	int size;
	fseek(fpbmp, 18L, SEEK_SET);
	fread(&width, sizeof(char), 4, fpbmp);
	fseek(fpbmp, 2L, SEEK_SET);
	fread(&size, sizeof(char), 4, fpbmp);
	printf("The Size of the bmp file is %ld.\n", size);
	fseek(fpbmp, 22L, SEEK_SET);
	fread(&height, sizeof(char), 4, fpbmp);
	printf("The Width of the bmp file is %ld.\n", width);
	printf("The Height of the bmp file is %ld.\n", height);
}

//=====================here is what I write=================

void RGBtoYUV()
{
	//convert
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			Y[j][i]= 0.299*r[j][i] + 0.587*g[j][i] + 0.114*b[j][i];
			U[j][i]= -0.147*r[j][i] - 0.289*g[j][i] + 0.436*b[j][i] ;
			V[j][i]= 0.615*r[j][i] - 0.515*g[j][i] - 0.100*b[j][i];
		}
	}
			
}
void YUVtoGray()//output YUV
{

	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				output_Y[j][i]=Y[j][i];
				output_U[j][i]=0;
				output_V[j][i]=0;
		}
	}
	
}


void YUVtoRGB()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				output_r[j][i] = output_Y[j][i] + 1.140*output_V[j][i];
				output_g[j][i] = output_Y[j][i] - 0.394*output_U[j][i] - 0.581*output_V[j][i];
				output_b[j][i] = output_Y[j][i] + 2.032*output_U[j][i];
		}
	}
		
}
void ChangeTheLuminance()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				output_Y[j][i]=255-Y[j][i];	
				output_U[j][i]=U[j][i];
				output_V[j][i]=V[j][i];
		}
	}
}

void YUVtoBinarization()//output YUV
{
		
	//OTSU althrithm '大津算法'
	unsigned char T=0,perfect_T=0;//Threshold
	double w0,w1,u0,u1,u;//paramenter
	long number_below_T = 0L,number_above_T = 0L;//paramenter
	double g[255]={0.0},max_g=0.0;//result;
	
	for(int num=0;num<=255;num++,T++)//after this there is a max g
	{
		/*Pre calcu*/
		long all_below_T=0L,all_above_T=0L;
		number_below_T=0L,number_above_T=0L;
		for (int j = 0; j<height; j++)
		{
			for (int i = 0; i<width; i++)
			{
					if(T>=Y[j][i]){
						number_below_T++;
						all_below_T+=Y[j][i];
					}
						
					else{
						number_above_T++;
						all_above_T+=Y[j][i];
					}
			}
		}
		if(number_below_T==0||number_above_T==0) continue;
		u0=(double)all_below_T/number_below_T;
		u1=(double)all_above_T/number_above_T;
		/*Step2*/
		w0=(double)number_below_T/(height*width);
		w1=(double)number_above_T/(height*width);
		u=w0*u0+w1*u1;
		g[num]=w0*pow(abs(u0-u),2.0)+w1*pow(abs(u1-u),2.0);
		/*Step 3 find max g*/
		if(g[num]>max_g){
			max_g=g[num];
			perfect_T=T;
		}
	}
	//'用算出来的perfect_T解出二值图'
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				if(Y[j][i]<perfect_T)
					output_Y[j][i]=0;
				else
					output_Y[j][i]=255;
				output_U[j][i]=0;
				output_V[j][i]=0;
		}
	}		
	printf("Binarization finished!\n");
}


void Dilation_Lining(char dowhat)
{
	int Struture_Element[elementline][elementline]={{0}};//假设是十字型
	int centre=(elementline-1)/2;
	/*自动生成十字形element*/
	for (int j = 0; j<elementline; j++)
	{
		for (int i = 0; i<elementline; i++)
		{
			if(j==centre||i==centre)
				Struture_Element[j][i]=1;
			else
				Struture_Element[j][i]=0;
		}
	}
	
		/*叠加原图*/
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				output_Y[j][i]=Y[j][i];
				output_U[j][i]=0;
				output_V[j][i]=0;
		}
	}
	
	/*膨胀一小圈*/
	for (int j = 0; j<height-elementline; j++)
	{
		for (int i = 0; i<width-elementline; i++)
		{
				int Istrue=0;
				for(int i1=0;i1<elementline;i1++)
					for(int j1=0;j1<elementline;j1++)
					{
						if(Struture_Element[i1][j1]==1&&Y[j+i1][i+j1]==0)
							Istrue=1;
					}
				/*大概，YUV中黑色是。。。Y=0 （Ｔ＾Ｔ）*/
				if(Istrue==1)
					output_Y[j+centre][i+centre]=dowhat=='d'?0:255;
				
		}
	}
	printf("Dilation successed!\n");
}
void Erosion()
{
	int Struture_Element[elementline][elementline]={{0}};//假设是十字型
	int centre=(elementline-1)/2;
	/*自动生成十字形element*/
	for (int j = 0; j<elementline; j++)
	{
		for (int i = 0; i<elementline; i++)
		{
			if(j==centre||i==centre)
				Struture_Element[j][i]=1;
			else
				Struture_Element[j][i]=0;
		}
	}
	
		/*叠加原图*/
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				output_Y[j][i]=255;
				output_U[j][i]=0;
				output_V[j][i]=0;
		}
	}
	
	/*缩小一小圈*/
	for (int j = 0; j<height-elementline; j++)
	{
		for (int i = 0; i<width-elementline; i++)
		{
				int Istrue=0;
				for(int i1=0;i1<elementline;i1++)
					for(int j1=0;j1<elementline;j1++)
					{
						if(Struture_Element[i1][j1]==0)
							continue;
						else if(Struture_Element[i1][j1]==1&&Y[j+i1][i+j1]!=0)//不同时黑色
							Istrue=1;//不全等
					}
				//printf("%d",Istrue); 
				/*大概，YUV中黑色是。。。Y=0 （Ｔ＾Ｔ）*/
				if(Istrue==0)
					output_Y[j+centre][i+centre]=0;//所有情况都成立！
				
		}
	}
	printf("Erosion successed!\n");
}
void Opening_Closing(char dowhat)
{
	if(dowhat=='o'){
		Erosion();
		BetweenTwoOperationProcessing();
		Dilation_Lining('d');
		printf("Opening finished, Now outer noises are removed!\n");
	}
	else if(dowhat=='c'){
		Dilation_Lining('d');
		BetweenTwoOperationProcessing();
		Erosion();
		printf("Closing finished, Now inner noises are removed!\n");
	}
	else
		printf("Opening_Closing fail!\n");
}
void BetweenTwoOperationProcessing()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				Y[j][i] = output_Y[j][i];
				U[j][i] = output_U[j][i];
				V[j][i] = output_V[j][i];
		}
	}
}
//============================3==============================
void VisibilityEnhancement()
{
	int maxLumnance = 0;
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				if(Y[j][i]>maxLumnance)
					maxLumnance=Y[j][i];
		}
	}
	
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				output_Y[j][i]=(unsigned char)255*(log10(Y[j][i]+1))/(log10(maxLumnance+1));	
				output_U[j][i]=0;
				output_V[j][i]=0;
		}
	}
	
}
void HistogramEqualization_r()
{
	//number 统计
	long n0[256]={0.0};
	/*0-255 直方图X3 RGB*/
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			n0[r[j][i]]++;
		}
	}

	/*0-1/255-...-254/255,1*/
	double p0[256]={0.0};
	for(int i=0;i<256;i++)
		p0[i]=n0[i]/(1.0*height*width);
		
	/*n[k]0-255:number*/
	
	/*Pn[n]*/
	double s0[256]={0.0};
	s0[0]=p0[0];
	for (int i=1;i<256;i++)
		s0[i]=s0[i-1]+p0[i];
//---------------------------------
	/*s[k]=求和Pn[0+..+k]*/
	/*找距离s[k]最近的像素点s[x]*/
	double minmize=1.0;
	int min_number[256]={0};
	double diff=0.0;
	for(int j = 0;j<256;j++)
	{
		//printf("!! == %f\n",s0[j]);
		diff=0.0;minmize=1.0;
		for(int i=0;i<256;i++)
		{
			double pix0=1.0*i/256;
			diff = fabs(s0[j]-pix0);
			//printf("%d == %f\n",i,diff);
			if(diff<minmize){
				minmize=diff;
				min_number[j]=i;//min_number[j] --> i ;
			}
				
		}
		//printf("%d == %d\n",j,min_number[j]);
	//	break;
	}
	//输出
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			output_r[j][i]=min_number[r[j][i]];
		}
	}
		
	/*统计临近的相同的像素点*/
	/*计算最终的而结果并画直方图*/
}void HistogramEqualization_g()
{
	//number 统计
	long n0[256]={0.0};
	/*0-255 直方图X3 RGB*/
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			n0[g[j][i]]++;
		}
	}

	/*0-1/255-...-254/255,1*/
	double p0[256]={0.0};
	for(int i=0;i<256;i++)
		p0[i]=n0[i]/(1.0*height*width);
		
	/*n[k]0-255:number*/
	
	/*Pn[n]*/
	double s0[256]={0.0};
	s0[0]=p0[0];
	for (int i=1;i<256;i++)
		s0[i]=s0[i-1]+p0[i];
//---------------------------------
	/*s[k]=求和Pn[0+..+k]*/
	/*找距离s[k]最近的像素点s[x]*/
	double minmize=1.0;
	int min_number[256]={0};
	double diff=0.0;
	for(int j = 0;j<256;j++)
	{
		//printf("!! == %f\n",s0[j]);
		diff=0.0;minmize=1.0;
		for(int i=0;i<256;i++)
		{
			double pix0=1.0*i/256;
			diff = fabs(s0[j]-pix0);
			//printf("%d == %f\n",i,diff);
			if(diff<minmize){
				minmize=diff;
				min_number[j]=i;//min_number[j] --> i ;
			}
				
		}
		//printf("%d == %d\n",j,min_number[j]);
	//	break;
	}
	//输出
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			output_g[j][i]=min_number[g[j][i]];
		}
	}
		
	/*统计临近的相同的像素点*/
	/*计算最终的而结果并画直方图*/
}void HistogramEqualization_b()
{
	//number 统计
	long n0[256]={0.0};
	/*0-255 直方图X3 RGB*/
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			n0[b[j][i]]++;
		}
	}

	/*0-1/255-...-254/255,1*/
	double p0[256]={0.0};
	for(int i=0;i<256;i++)
		p0[i]=n0[i]/(1.0*height*width);
		
	/*n[k]0-255:number*/
	
	/*Pn[n]*/
	double s0[256]={0.0};
	s0[0]=p0[0];
	for (int i=1;i<256;i++)
		s0[i]=s0[i-1]+p0[i];
//---------------------------------
	/*s[k]=求和Pn[0+..+k]*/
	/*找距离s[k]最近的像素点s[x]*/
	double minmize=1.0;
	int min_number[256]={0};
	double diff=0.0;
	for(int j = 0;j<256;j++)
	{
		//printf("!! == %f\n",s0[j]);
		diff=0.0;minmize=1.0;
		for(int i=0;i<256;i++)
		{
			double pix0=1.0*i/256;
			diff = fabs(s0[j]-pix0);
			//printf("%d == %f\n",i,diff);
			if(diff<minmize){
				minmize=diff;
				min_number[j]=i;//min_number[j] --> i ;
			}
				
		}
		//printf("%d == %d\n",j,min_number[j]);
	//	break;
	}
	//输出
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			output_b[j][i]=min_number[b[j][i]];
		}
	}
		
	/*统计临近的相同的像素点*/
	/*计算最终的而结果并画直方图*/
}

void HistogramEqualization_Y()
{
	//number 统计
	long n0[256]={0.0};
	/*0-255 直方图X3 RGB*/
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			n0[Y[j][i]]++;
		}
	}

	/*0-1/255-...-254/255,1*/
	double p0[256]={0.0};
	for(int i=0;i<256;i++)
		p0[i]=n0[i]/(1.0*height*width);
		
	/*n[k]0-255:number*/
	
	/*Pn[n]*/
	double s0[256]={0.0};
	s0[0]=p0[0];
	for (int i=1;i<256;i++)
		s0[i]=s0[i-1]+p0[i];
//---------------------------------
	/*s[k]=求和Pn[0+..+k]*/
	/*找距离s[k]最近的像素点s[x]*/
	double minmize=1.0;
	int min_number[256]={0};
	double diff=0.0;
	for(int j = 0;j<256;j++)
	{
		//printf("!! == %f\n",s0[j]);
		diff=0.0;minmize=1.0;
		for(int i=0;i<256;i++)
		{
			double pix0=1.0*i/256;
			diff = fabs(s0[j]-pix0);
			//printf("%d == %f\n",i,diff);
			if(diff<minmize){
				minmize=diff;
				min_number[j]=i;//min_number[j] --> i ;
			}
				
		}
		//printf("%d == %d\n",j,min_number[j]);
	//	break;
	}
	//输出
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			output_Y[j][i]=min_number[Y[j][i]];
		}
	}
		
	/*统计临近的相同的像素点*/
	/*计算最终的而结果并画直方图*/
}
//--------------------------4-----------------------------
void translation(int sizeX,int sizeY,int sizeX1,int sizeY1,int x,int y)
{
	int translatingMatrix[3][3] =  {{1,0,x},
									{0,1,y},
									{0,0,1}};
	//int x1 = translatingMatrix[0][0]*x0+translatingMatrix[0][1]*y0+translatingMatrix[0][2];
	//int y1 = translatingMatrix[1][0]*x0+translatingMatrix[1][1]*y0+translatingMatrix[1][2];
	if(sizeX>sizeX1||sizeY>sizeY1){
		printf("Be careful to mentain this: sizeX<=sizeX1 and sizeY<=sizeY1.\n");
		return;
	}
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
				output_r[y0][x0]= r[y0][x0];
				output_g[y0][x0]= g[y0][x0];
				output_b[y0][x0]= b[y0][x0];
			
				
		}
	}
	
	for (int y0 = sizeY; y0<sizeY1; y0++)
	{
		for (int x0 = sizeX; x0<sizeX1; x0++)
		{
			int x1 = translatingMatrix[0][0]*x0+translatingMatrix[0][1]*y0+translatingMatrix[0][2];
			int y1 = translatingMatrix[1][0]*x0+translatingMatrix[1][1]*y0+translatingMatrix[1][2];
			
			output_r[y1][x1]= r[y0][x0];
			output_g[y1][x1]= g[y0][x0];
			output_b[y1][x1]= b[y0][x0];
		}
	}
	for (int y0 = sizeY+y; y0<sizeY1+y; y0++)
	{
		for (int x0 = sizeX+x; x0<sizeX1+x; x0++)
		{
			int x1 = translatingMatrix[0][0]*x0+translatingMatrix[0][1]*y0-translatingMatrix[0][2];
			int y1 = translatingMatrix[1][0]*x0+translatingMatrix[1][1]*y0-translatingMatrix[1][2];
			
			output_r[y1][x1]= r[y0][x0];
			output_g[y1][x1]= g[y0][x0];
			output_b[y1][x1]= b[y0][x0];
		}
	}
	
}
void translation1(int sizeX,int sizeY,int sizeX1,int sizeY1,int x,int y)
{
	int translatingMatrix[3][3] =  {{1,0,x},
									{0,1,y},
									{0,0,1}};
	//int x1 = translatingMatrix[0][0]*x0+translatingMatrix[0][1]*y0+translatingMatrix[0][2];
	//int y1 = translatingMatrix[1][0]*x0+translatingMatrix[1][1]*y0+translatingMatrix[1][2];
	if(sizeX>sizeX1||sizeY>sizeY1){
		printf("Be careful to mentain this: sizeX<=sizeX1 and sizeY<=sizeY1.\n");
		return;
	}
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
				output_r[y0][x0]= 255;
				output_g[y0][x0]= 255;
				output_b[y0][x0]= 255;				
		}
	}
	
	for (int y0 = sizeY; y0<sizeY1; y0++)
	{
		for (int x0 = sizeX; x0<sizeX1; x0++)
		{
			int x1 = translatingMatrix[0][0]*x0+translatingMatrix[0][1]*y0+translatingMatrix[0][2];
			int y1 = translatingMatrix[1][0]*x0+translatingMatrix[1][1]*y0+translatingMatrix[1][2];
			
			output_r[y1][x1]= r[y0][x0];
			output_g[y1][x1]= g[y0][x0];
			output_b[y1][x1]= b[y0][x0];
		}
	}	
}
void mirror(char character)
{
	if(character=='x')
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
				output_r[y0][x0]= r[y0][width-x0];
				output_g[y0][x0]= g[y0][width-x0];
				output_b[y0][x0]= b[y0][width-x0];	
		}
	}
	else if(character=='y')
		for (int y0 = 0; y0<height; y0++)
		{
			for (int x0 = 0; x0<width; x0++)
			{
					output_r[y0][x0]= r[height-y0][x0];
					output_g[y0][x0]= g[height-y0][x0];
					output_b[y0][x0]= b[height-y0][x0];	
			}
		}
	else
		printf("something wrong with your input.");
}
void Scaling(float ratioX,float ratioY)
{
	if(ratioX<=0||ratioY<=0){
		printf("Scaling paramenter wrong!\n");
		return;
	}
	if(ratioX<=1&&ratioY<=1){
		for (int y0 = 0; y0<height; y0++)
		{
			for (int x0 = 0; x0<width; x0++)
			{
				int x1 = ratioX*x0;
				int y1 = ratioY*y0;
				
				output_r[y1][x1]= r[y0][x0];
				output_g[y1][x1]= g[y0][x0];
				output_b[y1][x1]= b[y0][x0];
				
			}
		}
	}
	else if(ratioX>1||ratioY>1){
		
		for (int y0 = 0; y0<height; y0++)
		{
			for (int x0 = 0; x0<width; x0++)
			{
				int x1 = ratioX*x0;
				int y1 = ratioY*y0;
				
				output_r[y1][x1]= r[y0][x0];
				output_g[y1][x1]= g[y0][x0];
				output_b[y1][x1]= b[y0][x0];
				records[y1][x1]=1;
			}
		}
		
		for (int y0 = 0; y0<height*ratioY; y0++)
		{
			for (int x0 = 0; x0<width*ratioX; x0++)
			{
				if(records[y0][x0]==0){
					if(ratioX>1&&ratioY>1){
						output_r[y0][x0]= r[(int)(y0/ratioY)][(int)(x0/ratioX)];
						output_g[y0][x0]= g[(int)(y0/ratioY)][(int)(x0/ratioX)];
						output_b[y0][x0]= b[(int)(y0/ratioY)][(int)(x0/ratioX)];
					}
					else if(ratioX<=1&&ratioY>1){
						output_r[y0][x0]= r[(int)(y0/ratioY)][x0];
						output_g[y0][x0]= g[(int)(y0/ratioY)][x0];
						output_b[y0][x0]= b[(int)(y0/ratioY)][x0];
						
					}
					else if(ratioX>1&&ratioY<=1){
						output_r[y0][x0]= r[y0][(int)(x0/ratioX)];
						output_g[y0][x0]= g[y0][(int)(x0/ratioX)];
						output_b[y0][x0]= b[y0][(int)(x0/ratioX)];
						
					}
				}
				
			}
		}
	}
	
}
double Gaussian(int radial,double lamda)
{
	return exp(-pow(radial,2)/(2*pow(lamda,2)));	
}

void RBF_Gaussian(unsigned char temple[],int x,int y,int paramenter)//3
{
	/*Paramenter is used to control how many points are included included in this process*/
	//L2 fomular
	double lamda = 0.4;//paramenter_X
	double w0 = 0.0;
	for (int y0 = y-paramenter; y0<=y+paramenter; y0++)
		{
			for (int x0 = x-paramenter; x0<=x+paramenter; x0++)
			{
					if(x0==x&&y0==y) continue;
					double radial = abs(x-x0)+abs(y-y0);//sqrt(pow((x-x0),2)+pow((y-y0),2));
					w0 += Gaussian(radial,lamda);
			}
		}
	double pixelout[3] = {0.0};
	for (int y0 = y-paramenter; y0<=y+paramenter; y0++)
		{
			for (int x0 = x-paramenter; x0<=x+paramenter; x0++)
			{
					if(x0==x&&y0==y) 
						continue;
					double radial = abs(x-x0)+abs(y-y0);//sqrt(pow((x-x0),2)+pow((y-y0),2));
					double w1 = Gaussian(radial,lamda)/w0;
					pixelout[0] += w1*output_r[y0][x0];
					pixelout[1] += w1*output_g[y0][x0];
					pixelout[2] += w1*output_b[y0][x0];
			}
		}
	for(int i=0;i<3;i++) temple[i] = (unsigned char)pixelout[i];//output!
	//printf(" %d\n",temple[0]);
}


void Rotation(double angle)
{
	angle = (angle-90)/180*3.1415;
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
			//int x1 = cos(angle)*x0-sin(angle)*y0;
			//int y1 = sin(angle)*x0+cos(angle)*y0;
			int x1 = sin(angle)*(width/2-x0)+cos(angle)*(height/2-y0);
			int y1 = sin(angle)*(height/2-y0)-cos(angle)*(width/2-x0);
			y1+=height/2;x1+=width/2;
			records[y1][x1] = 1;
			output_r[y1][x1]= r[y0][x0];
			output_g[y1][x1]= g[y0][x0];
			output_b[y1][x1]= b[y0][x0];
		}
	}
	
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
			if(records[y0][x0] != 1)
			{		
				/*指针传递三个参数值*/
				unsigned char temple[3];
				RBF_Gaussian(temple,x0,y0,3);
				output_r[y0][x0]=temple[0];
				output_g[y0][x0]=temple[1];
				output_b[y0][x0]=temple[2];
			}
		}
	}
	
}
void Shearing(float dx,float dy)
{
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
			int x1=0; 
			int y1=0;
			if(dx==0&&dy!=0){
				x1=x0;
				y1=y0+dy*x0;
			}
			else if(dx!=0&&dy==0){
				x1 = x0+dx*y0;
				y1 = y0;
			}
			else if(dx!=0&&dy!=0){
				x1 = x0+dx*y0;
				y1 = y0+dy*x0;
			}
			else
				printf("Shearing paramenter errors.\n");
			
			output_r[y1][x1]= r[y0][x0];
			output_g[y1][x1]= g[y0][x0];
			output_b[y1][x1]= b[y0][x0];
		}
	}
	
}
//------------------------------5--------------------------------
void MeanFilter(unsigned char temple[],int x,int y,int paramenter)//给出该点的准确值=1
{
		double wr=0.0,wg=0.0,wb=0.0;
		int count=0;			
		for (int y0 = y-paramenter; y0<=y+paramenter; y0++)
		{
			for (int x0 = x-paramenter; x0<=x+paramenter; x0++)
			{
					if(x0<0||x0>width-1||y0<0||y0>height-1)
						continue;
					wr+=r[y0][x0];
					wg+=g[y0][x0];
					wb+=b[y0][x0];
					count++;
			}
		}
		wr/=count;wg/=count;wb/=count;
		temple[0]=(unsigned char)wr;
		temple[1]=(unsigned char)wg;
		temple[2]=(unsigned char)wb;
}

void meanFiltering(int paramenter)
{
	
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
				unsigned char temple[3];
				MeanFilter(temple,x0,y0,paramenter);
				output_r[y0][x0]=temple[0];
				output_g[y0][x0]=temple[1];
				output_b[y0][x0]=temple[2];
		}
	}
	
	printf("Mean Filtering Success!\n");
	
}
void Meanfilter_Y(int paramenter)
{
	for (int y = 0; y<height; y++)
	{
		for (int x = 0; x<width; x++)
		{
			double w0=0.0;
			int count=0;			
			for (int y0 = y-paramenter; y0<=y+paramenter; y0++)
			{
				for (int x0 = x-paramenter; x0<=x+paramenter; x0++)
				{
					if(x0<0||x0>width-1||y0<0||y0>height-1)
						continue;
					w0+=Y[y0][x0];
					count++;
				}
			}
			output_Y[y][x]=w0/count;
		}
	}	
	
	
}

unsigned char Laplacian_filter_Y(int x,int y)
{
	int paramenter=1;
	int w0=0;		
		for (int y0 = y-paramenter; y0<=y+paramenter; y0++)
		{
			for (int x0 = x-paramenter; x0<=x+paramenter; x0++)
			{
					if(x0<0||x0>width-1||y0<0||y0>height-1)
						continue;
					if(x==x0&&y==y0){
						w0-=4*Y[y0][x0];
					}
					else if(abs(x+y-x0-y0)==1){
						w0+=Y[y0][x0];
					}
			}
		}
	return (unsigned char)w0;
}

void LaplacianFilter_Y(double paramenter)
{
	
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{
			output_Y[y0][x0]=0;
			output_Y[y0][x0]+=abs(paramenter*Laplacian_filter_Y(x0,y0));//nomal
			
			output_Y[y0][x0]+=Y[y0][x0];
		}
	}
	
	printf("Laplacian Filtering Success!\n");
	
}

//--------------------------6--------------------------------
void RBF_Bilateral(unsigned char temple[],int x,int y,int paramenter,double lamda_range,double lamda_color)
{
	double w0_r = 0.0;
	double w0_g = 0.0;
	double w0_b = 0.0;
	for (int y0 = y-paramenter; y0<=y+paramenter; y0++)
		{
			for (int x0 = x-paramenter; x0<=x+paramenter; x0++)
			{
					if(x0==x&&y0==y) continue;
					double radial = abs(x-x0)+abs(y-y0);//sqrt(pow((x-x0),2)+pow((y-y0),2));
					double color_r = abs(r[y][x]-r[y0][x0]);
					double color_g = abs(g[y][x]-g[y0][x0]);
					double color_b = abs(b[y][x]-b[y0][x0]);
					w0_r += Gaussian(radial,lamda_range)*Gaussian(color_r,lamda_color);
					w0_g += Gaussian(radial,lamda_range)*Gaussian(color_g,lamda_color);
					w0_b += Gaussian(radial,lamda_range)*Gaussian(color_b,lamda_color);
			}
		}
	double pixelout[3] = {0.0};
	for (int y0 = y-paramenter; y0<=y+paramenter; y0++)
		{
			for (int x0 = x-paramenter; x0<=x+paramenter; x0++)
			{
					if(x0==x&&y0==y) continue;
					double radial = abs(x-x0)+abs(y-y0);
					double color_r = abs(r[y][x]-r[y0][x0]);
					double color_g = abs(g[y][x]-g[y0][x0]);
					double color_b = abs(b[y][x]-b[y0][x0]);
					double w1_r = Gaussian(radial,lamda_range)*Gaussian(color_r,lamda_color)/w0_r;
					double w1_g = Gaussian(radial,lamda_range)*Gaussian(color_g,lamda_color)/w0_g;
					double w1_b = Gaussian(radial,lamda_range)*Gaussian(color_b,lamda_color)/w0_b;
					
					pixelout[0] += w1_r*r[y0][x0];
					pixelout[1] += w1_g*g[y0][x0];
					pixelout[2] += w1_b*b[y0][x0];
			}
		}
		
	for(int i=0;i<3;i++) temple[i] = (unsigned char)pixelout[i];//output!
	
}
void BilateralFilter(int paramenter,double lamda_range,double lamda_color)
{
	clock_t start,finish;
	start=clock();
	for (int y0 = 0; y0<height; y0++)
	{
		for (int x0 = 0; x0<width; x0++)
		{		
				/*指针传递三个参数值*/
				unsigned char temple[3];
				RBF_Bilateral(temple,x0,y0,paramenter,lamda_range,lamda_color);
				output_r[y0][x0]=temple[0];
				output_g[y0][x0]=temple[1];
				output_b[y0][x0]=temple[2];
				//printf("%d\n",temple[0]-r[y0][x0]);
		}
	}
	finish=clock();
	printf("Bilateral Filter finished.\n");
	int duration=(int)(finish-start)/CLOCKS_PER_SEC;
	printf("\aIt costs %02d:%02d \n",duration/60,duration%60);
}






//===========================End=================================

//===========================改写文件流==========================
void Initialization()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
				output_Y[j][i]=255;
				output_U[j][i]=0;
				output_V[j][i]=0;
				records[j][i]=0;
		}
	}
		
}
//改写文件流
FILE *openfile(char filename[])
{
	FILE *fpbmp;
	fpbmp = fopen(filename, "rb");
	if (fpbmp == NULL)
	{
		printf("Open bmp failed!!!\n");
		exit(1);
	}
	
	bmpFileTest(fpbmp);
	bmpHeaderPartLength(fpbmp);
	BmpWidthHeight(fpbmp);
	
	fseek(fpbmp, 0L, SEEK_SET);
	return fpbmp;
}

FILE *writefile(char filename[])
{
	FILE *fpout;
	fpout = fopen(filename, "wb+");
	if (fpout == NULL)
	{
		printf("Open out.bmp failed!!!\n");
		exit(1);
	}
	fseek(fpout, 0L, SEEK_SET);
	return fpout;
	
}

void addHeadertofile(FILE *input,FILE *output)
{
		unsigned char *fp_temp;
		
		fseek(input, 0L, SEEK_SET);
		fseek(output, 0L, SEEK_SET);
		
		fp_temp = malloc(OffSet);
		fread(fp_temp, 1, OffSet, input);//读取头文件
		fp_temp[18]=(int)width;
		fp_temp[22]=(int)height;
		fp_temp[2]=(int)(OffSet+height*width*3);
		fp_temp[34]=(int)(height*((24 * width/8 + 3)/4*4));
		//printf("%#X %#X %d",fp_temp[2],fp_temp[4],fp_temp[0]);
		fwrite(fp_temp, 1, OffSet, output);//输出头文件
		

}
//=====================改写文件流（完）================ 

//=====================TEST===========================
void testInputYUV()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			printf("%d ",Y[j][i]);
		}
	}
}
void testOutputYUV()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			printf("%d ",output_Y[j][i]);
		}
	}
}
void testRGB()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			printf("%d %d %d \n",r[j][i],g[j][i],b[j][i]);
		}
	}	
}
void NoneProcess()
{
	for (int j = 0; j<height; j++)
	{
		for (int i = 0; i<width; i++)
		{
			output_r[j][i]= r[j][i];
			output_g[j][i]= g[j][i];
			output_b[j][i]= b[j][i];
		}
	}	
}
//=====================TEST END=================
