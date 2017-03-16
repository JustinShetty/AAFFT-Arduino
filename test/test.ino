/*
REAL FFT tests.
atmega1280 @ 16 MHz with 32k external sram
910 milliseconds for 2048 point floating point transform (100 ms/256 points)
// relocated heap to external SRAM via linker command
// -Wl,--defsym=__heap_start=0x801100,--defsym=__heap_end=0x807fff
// also must include floating point printf and lib:
// -Wl,-u,vfprintf -lprintf_flt -lm

*/

#include <avr/io.h>
#include <avr/interrupt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <avr/pgmspace.h>
#define F_CPU 11059200UL
#include <util/delay.h>      


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define NP 16

///////////////////////////////////////////////////////
// Global Variables

float *data;	//allocated at runtime in main
int ndata=NP;	//size of transform buffer

// initialization sequence in .init3

void init_extSRAM(void) __attribute__((naked,section(".init3"))) ;

void init_extSRAM(void) {
//Enable external SRAM: set MCUCR, XMCRA and XMCRB for
// External SRAM page configuration: 0000h - 7FFFh
// Lower page wait state(s): None
// Upper page wait state(s): None
	MCUCR=0x80;
	XMCRA=0x00;    
	XMCRB &=~7 ;  //Make bits 0-2 of XMCRB low to enable all of external SRAM
}

volatile unsigned short ticks=0;	//0.01 second clock ticks updated in ISR

//  delay for time_ms milliseconds by looping
// _delay_ms() has max delay about 13 ms @ 20 MHz

void delay_ms( unsigned short time_ms )
{
	unsigned int i;
	for ( i = 0; i < time_ms; i++ )  _delay_ms( 1 );
}

// global COM=0 or 1 specifies UART0 or UART1

int COM=1;

//#include "uputget.c"

//int uart_putchar(char c, FILE *stream) {
/*
* Send character c to the LCD for use with printf
*/
//   uputchar(c);
//   return 0;
// }

//setup outpuf file for printf

//FILE uart_str = FDEV_SETUP_STREAM(uart_putchar, NULL, _FDEV_SETUP_WRITE);

/*
* Timer0 overflow interrupt handler.  
* provides global ticks of 10 ms
*/

ISR(TIMER0_OVF_vect)
{
	TCNT0=256-108;
	ticks++; 		
}

void timer_init(void) {

// set up Timer0 to provide 100 tick/sec

//	TCCR0  = (7<<CS00);  //clock/1024
//	TIMSK |= (1<<TOIE0); // enable int on overflow
//	TCNT1 = 256-108; //exact divisor for 11.059 MHz
}

void get_data(void)
{
	int j;
	for (j=0;j<NP/2;j++) {
	data[j]= -1.0;
	data[j+NP/2]= 1.0;
	}
}
void four1(int nn, int isign)
{
	int n,mmax,m,j,istep,i;
	/*double*/
	float wtemp,wr,wpr,wpi,wi,theta;

	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.2831853/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void realft(int n, int isign)
{

	int i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	/*double*/
	float wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592654/(float) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(n>>1,-1);
	}
}

void printDouble( double val, unsigned int precision){
// prints val with number of decimal places determine by precision
// NOTE: precision is 1 followed by the number of zeros for the desired number of decimial places
// example: printDouble( 3.1415, 100); // prints 3.14 (two decimal places)

   Serial.print (int(val));  //prints the int part
   Serial.print("."); // print the decimal point
   unsigned int frac;
   if(val >= 0)
     frac = (val - int(val)) * precision;
   else
      frac = (int(val)- val ) * precision;
   int frac1 = frac;
   while( frac1 /= 10 )
       precision /= 10;
   precision /= 10;
   while(  precision /= 10)
       Serial.print("0");

   Serial.print(frac,DEC) ;
}

#include <math.h>

short FFT(short int dir,int m,float *x,float *y)
{
   int n,i,i1,j,k,i2,l,l1,l2;
   float c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }

   return(0);
}

////////////////////////////////////////////////////////////////////////
int main ()
{	
  Serial.begin(9600);
	int i;
	DDRD |= 1;  	//set PORTD.0 = output
	PORTD &= ~1;  	//led on
	delay_ms(200);
	PORTD |= 1;	//led off
		
	//uart_init();            
	//stdout = &uart_str;         // associate stream with stdout
	//timer_init();
	sei();

	Serial.print("Real FFT 0.1\n\r");   // announce string

	if ( (data = (float *)malloc( (ndata+1)*sizeof(float)) ) == NULL )
		Serial.print("Unable to allocate data array\r\n");
	else
		Serial.print("Successfully allocated data array\r\n");
		
	get_data();

	for(i=0; i<16; i+=2){
		//printf( "%7.2f %7.2f \r\n",(double) data[i],(double) data[i+1]);
		printDouble(data[i],100000);
		//Serial.print(data[i]);
		Serial.print(" ");
    printDouble(data[i+1],100000);
		//Serial.println(data[i+1]);
    Serial.println(" ");
	}
	ticks=0;
	short FFT(1,4,data,data);
	realft(NP,1); //transform
	i=ticks;
	//printf(" realfft in %dx10 ms \r\n",i);
	Serial.print(NP);
	Serial.println(" data points \r\n");
	//printf(" %d data points \r\n",NP);
	for(i=0; i<16; i+=2){
		//printf( "%7.2f %7.2f \r\n",(double) data[i],(double) data[i+1]);
    printDouble(data[i],100000);
		//Serial.print(data[i]);
		Serial.print(" ");
    printDouble(data[i+1],100000);
		//Serial.println(data[i+1]);
    Serial.println(" ");
	}
	while(1);
}
