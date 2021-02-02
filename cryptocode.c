
#include<stdio.h>
#include<gmp.h>
#include <string.h>
struct Elliptic_Curve   //structure containing parameters of curve
{
mpz_t a;
mpz_t b;
mpz_t p;
};
struct Point       //structure containig Point values
{
mpz_t x;
mpz_t y;
};
struct Elliptic_Curve EC;

//Set parameters a , b and c in ... y^2 = ( x^3 + ax + b ) mod p
void Select_EC()
{
printf("\nEnter Elliptic Curve Parameters i.e. a, b and p: \n");
gmp_scanf("%Zd",&EC.a);
gmp_scanf("%Zd",&EC.b);
gmp_scanf("%Zd",&EC.p);
}

//performing R = P + P 
void Point_Doubling(struct Point P, struct Point *R)
{
mpz_t slope,temp;        
mpz_init(temp);       //slope = 0
mpz_init(slope);          //temp = 0


//Equation of slope of tangent : x = 3*x^2 + a / 2*y

if(mpz_cmp_ui(P.y,0)!=0)
{
mpz_mul_ui(temp,P.y,2) ;     //temp = 2 * y
mpz_invert(temp,temp,EC.p);  // temp = 1 / temp
mpz_mul(slope,P.x,P.x);     //slope = x^2 
mpz_mul_ui(slope,slope,3);   // slope = 3 * x^2

mpz_add(slope,slope,EC.a);   //slope = slope + a
mpz_mul(slope,slope,temp);   // slope = slope / temp
mpz_mod(slope,slope,EC.p);  // slope mod p

//Using diophantus method to plot points  i.e x1 + x2 + x = m^2 
//Here x1 = x2 ... x = m^2 - 2*x1

mpz_mul(R->x,slope,slope); // R->x = m^2
mpz_sub(R->x,R->x,P.x);     //R->x = m^2 - x
mpz_sub(R->x,R->x,P.x);   //R->x = m^2 - 2x
mpz_mod(R->x,R->x,EC.p);  //R->x mod p
//R->x plotted

//plotting R->y      .... using eqn : y1 = m^2 ( x2-x1 ) - y2
mpz_sub(temp,P.x,R->x);    //temp = x2 - x1
mpz_mul(R->y,slope,temp);  //y1 = m^2 * temp
mpz_sub(R->y,R->y,P.y);  //y1 = y1 - y2
mpz_mod(R->y,R->y,EC.p);  // y1 mod P
//R-y plotted
}


// if P.y == 0 , then P.x is on X-axis => tangent will meet at point O i.e infinity(0,0)
else
{
mpz_set_ui(R->x,0);
mpz_set_ui(R->y,0);
}
}


//R = P + Q
void Point_Addition(struct Point P,struct Point Q, struct Point*R)
{ 
mpz_mod(P.x,P.x,EC.p);
mpz_mod(P.y,P.y,EC.p);
mpz_mod(Q.x,Q.x,EC.p);
mpz_mod(Q.y,Q.y,EC.p);
mpz_t temp,slope;
mpz_init(temp);
mpz_init_set_ui(slope,0);


if(mpz_cmp_ui(P.x,0)==0 && mpz_cmp_ui(P.y,0)==0)  //if P = (0,0) , return Q
{ 
mpz_set(R->x,Q.x); 
mpz_set(R->y,Q.y); 
return ;
}
if(mpz_cmp_ui(Q.x,0)==0 && mpz_cmp_ui(Q.y,0)==0)   //if Q = (0,0), return P
{ 
mpz_set(R->x,P.x); 
mpz_set(R->y,P.y);
return ;
}

if(mpz_cmp(P.x,Q.x)==0 && mpz_cmp(P.y,Q.y)==0)   //if P==Q ... perform R = 2P
{
Point_Doubling(P,R);  // as point P = Q...perform 2P
return ;
}
else
{
mpz_sub(temp,P.x,Q.x);       // temp = x2 -x1
mpz_mod(temp,temp,EC.p);       //temp = temp mod p
mpz_invert(temp,temp,EC.p);   //temp = 1/temp
mpz_sub(slope,P.y,Q.y);       //slope = y2 -y1
mpz_mul(slope,slope,temp);    // slope = slope * temp
mpz_mod(slope,slope,EC.p);    //slope = slope mod p

//Plotting R->x
//Using Diophantus method of plotting points from known points i.e y = m^2 - x1 - x2
mpz_mul(R->x,slope,slope);  // m^2
mpz_sub(R->x,R->x,P.x);     // m^2 - x1
mpz_sub(R->x,R->x,Q.x);     // m^2 - x1 - x2  
mpz_mod(R->x,R->x,EC.p);    //R->x = R->x mod p.......R->x plotted 



//Now plotting R->y  by using  R->y = m(x-x1) - y

mpz_sub(temp,P.x,R->x);     //temp = P.x - R.x   
mpz_mul(R->y,slope,temp);  //R->y = m * temp
mpz_sub(R->y,R->y,P.y);    //R->y = R->y - P.y 
mpz_mod(R->y,R->y,EC.p);   //R->y mod p................R->y plotted

return ;
}
}

// R = mP 
void Scalar_Multiplication(struct Point P,struct Point *R,mpz_t m)
{
struct Point Q,T;   
mpz_init(Q.x);       
mpz_init(Q.y);
mpz_init(T.x); 
mpz_init(T.y);

long no_of_bits,loop;           
no_of_bits = mpz_sizeinbase(m,2);   //no. of digits in binary form of m
mpz_set_ui(R->x,0);
mpz_set_ui(R->y,0);

if(mpz_cmp_ui(m,0)==0)   //if slope = 0 .. O ponit i.e (0,0) will be returned
return;

mpz_set(Q.x,P.x);
mpz_set(Q.y,P.y);

if(mpz_tstbit(m,0)==1)  //if 0th bit of m == 1, 
{
mpz_set(R->x,P.x);
mpz_set(R->y,P.y);
}

// from the implementation of ECC curve, 
// loop will run (no of bits in m - 1 times)
// if bit encountered == 0 , perform Point_doubling i.e R -> 2P
// if bit encountered == 1 , perform Point_doubling then Point_Addition i.r R-> 2P + P
for(loop=1;loop<no_of_bits;loop++)      
{
mpz_set_ui(T.x,0);
mpz_set_ui(T.y,0);
Point_Doubling(Q,&T);
mpz_set(Q.x,T.x);
mpz_set(Q.y,T.y);
mpz_set(T.x,R->x);
mpz_set(T.y,R->y);

if(mpz_tstbit(m,loop))
Point_Addition(T,Q,R);

}
}

void hexa(mpz_t m){
    char p[1000];
    mpz_get_str(p,16,m);
    printf("%s\n",p);
    mpz_clear(m);
}


void main()
{
int choice;
mpz_init(EC.a);   //initializing EC.a
mpz_init(EC.b);   //initailizing EC.b
mpz_init(EC.p);   //initailizing EC.p i.e the GF(p)
Select_EC();         //setting values of a,b , p
printf("Enter Choice of Operation \n");
printf("\nPress '1' For Point Addition Operation \n");
printf("\nPress '2' For Scalar Multiplication Operation\n");
scanf("%d",&choice);
struct Point P,R;         
mpz_init(P.x);        //initializing P.x
mpz_init(P.y);            //initailizing P.y
mpz_init_set_ui(R.x,0);  //initializing and setting R.x = 0
mpz_init_set_ui(R.y,0);  //initializing and setting R.y = 0

printf("\nEnter Point P(x,y)\n");
gmp_scanf("%Zd",&P.x);  //user-input value
gmp_scanf("%Zd",&P.y);   //user input-value


if(choice==1)
{
struct Point Q;
printf("\nEnter Point Q(x,y)\n");
mpz_init(Q.x);   //initailizing Q.x
mpz_init(Q.y);    //initailizing Q.y


gmp_scanf("%Zd",&Q.x);   //user-input value
gmp_scanf("%Zd",&Q.y);  //user-input value
 

Point_Addition(P,Q,&R);
}
else
{
printf("\nEnter m to Find mP\n");
mpz_t m;
mpz_init(m);
gmp_scanf("%Zd",&m);
Scalar_Multiplication(P,&R,m);
}
gmp_printf("\nResultant Point is %Zd,%Zd\n\n",R.x,R.y);
printf("In the hexadecimal form:\n");
hexa(R.x);
hexa(R.y);
}
