

// Global variables.  These are used in IRFunctions.h and must be declared
//    before IRFunctions.h is included.

// Market data.
double *par, sigma;

// Lattice data. The usual stuff for the binomial lattice.
double **d;  // The state-dependent single-period discount factors.
double **V;  // State-dependent lattice values of future cash flow.

// Funky Bond data.
int    ***I;
double ***FunkyBond;
double    coupon;

#include "Functions.h"
#incle "IRFunctions.h"


// Functions found below the main program:
double ValueFunkyBond (int, int, int);
double ValueBond(int, int, int);
void   AllocateMemory ();



int main () {

   double left, right, middle, value, P1,P2,
         duration, faircoupon, fbond; 
   
   // Allocate all array space.
   AllocateMemory ();

   // Get the user specified par curve.
   MakeParCurve ();

   // Get the volatility parameter and convert to a decimal.
   sigma = GetDouble ("What is the short-rate volatility in percent?... ");
   sigma /= 100.0; // Now as a decimal.

   // Calibrate the Salomon binomial lattice to explain the par yield curve in par[].
   Calibrate ();

   printf ("\n");
   
   //finds fair coupon, initializes starting values
   left = 0.0;
   right = 100.0;
   middle = (left+right)/2.0; 
   coupon = middle; 
   //gets semiannual coupon as a decimal 
   coupon/= 200.0; 
   //calculates the difference in funkybond value and par 
   value = ValueFunkyBond(0,0, 20)-100.0;
   

   //if diff between funkybond and par is greater than 1 cent
   //uses interval bisection to find fair coupon rate
   while (abs(value)> .01){
      if (value >.01){
         right=middle;
      }

      if (value < -0.01){
         left = middle; 
      }
      
      middle = (left+right)/2.0;
      coupon = middle/200.0;  
      value = ValueFunkyBond(0,0, 20)-100.0;
   }
   
   //prints out fair coupon rate 
   printf("The fair coupon rate in percent is %8.3f\n", coupon*200); 
   //prints out value of bond with no deferred or doubling up on payments w/ fair coupon rate
   printf ("The Regular Bond's value is %8.3f\n", ValueBond (0,0,20));

   fbond = ValueFunkyBond(0,0,20);
   //prints out value of funky bond using fair coupon rate 
   printf ("The Funky Bond's value is %8.3f\n", fbond);
   printf("The combined value of the options to double up and/or defer principal payments %8.3f\n",
            ValueBond (0,0,20) - ValueFunkyBond (0,0,20) ); 
   printf("\n");
   

   //shifts par curve 
   AllocateMemory ();
   MakeParCurve ();
   Calibrate ();
   printf ("\n");
   //gets value of FunkyBond when par curve is shifted up by 25 basis points
   P1 = ValueFunkyBond(0,0,20); 

   //shifts par curve 
   AllocateMemory ();
   MakeParCurve ();
   Calibrate ();
   printf ("\n");
   //gets value of FunkyBond when par curve is shifted down by 25 basis points
   P2 = ValueFunkyBond(0,0,20);
   
   //calculates duration of FunkyBond 
   duration = -1.0/fbond * (P1-P2)/.005;
   printf("The Funky Bond's duration is %8.3f \n", duration); 



   // Pause so the execution window does not close.
   Pause ();

}

//values bond if principal payments are constant at $5 mil after period 20 
double ValueBond(int n, int i, int j){
   double c1; 
   //when first called, sets all indicator variables to be zero
   if (n == 0) {
      for (int n0 = 0; n0 <= 40; n0++) {
         for (int i0 = -n0; i0 <= n0; i0 += 2) {
            for (int j0 = 0; j0 <= 20; j0 ++) {
                  I[n0][i0][j0]= 0;
            }
         }
      }
   }

   //if no remaining principal -- value of bond is zero 
   if (j == 0){
      FunkyBond[n][i][j] = 0.0; 
      I[n][i][j] = 1; 
      }

   //if calculations have not been done yet at (n, i, j)
   //calculates future cash flows 
   if (I[n][i][j] == 0){

      //when in periods 1-19, there is only coupon payments and no principal payment
      if (n < 20){
         c1 = coupon*5*j;
         FunkyBond[n][i][j] = d[n][i] * ( 0.5 * (c1 + ValueBond(n+1, i+1,j)) + 
                           0.5* (c1 + ValueBond(n+1, i-1, j))); 
      }
      //calculates future cash flows of coupon and $5 principal payments from year 
      else{
         c1 = coupon*5*j + 5; 
         FunkyBond[n][i][j] = d[n][i] * ( 0.5 * (c1 + ValueBond(n+1, i+1,j-1)) + 
                           0.5* (c1 + ValueBond(n+1, i-1, j-1))); 
      }  
   }

   I[n][i][j] = 1; 
   return FunkyBond[n][i][j]; 

}

////////////////////////////////////////////////////////////////////////////////
// This function values a 20 year Funky Bond. Coupons are paid
//   semiannually. The bond is valued recursively.
// (n,i) is the node on the interest rate lattice and 5*j is the amount
//   of the bond outstanding in millions at period n AFTER any period n
//   principal payment.
////////////////////////////////////////////////////////////////////////////////
double ValueFunkyBond (int n, int i, int j) {

   double min, c1,c2, value;  
   min = 100000.0;
   //when first called, sets all indicator variables to be zero
   if (n == 0) {
      for (int n0 = 0; n0 <= 40; n0++) {
         for (int i0 = -n0; i0 <= n0; i0 += 2) {
            for (int j0 = 0; j0 <= 20; j0 ++) {
                  I[n0][i0][j0]= 0;
            }
         }
      }
   }

   //sets boundary condition to stop recursion
   //when there is no more remaining principal 
   if (j == 0 ) {
      FunkyBond[n][i][j] = 0.0; 
      I[n][i][j]= 1;}


   //if calculations have not been done at (n,i,j) yet
   if (I[n][i][j] == 0 ){
     

     //when we are in periods 1 -19 calculates future cash flow
     //only pays coupon * remaining principal in current period
     if (n < 20) {
      c1= coupon*5*j; 
      FunkyBond[n][i][j] = d[n][i] * ( 0.5 * (c1 + ValueFunkyBond(n+1, i+1,j)) + 
                           0.5* (c1 + ValueFunkyBond(n+1, i-1, j))); 
     }

      //calculates future cash flow for periods 20-40
      else{
         //cases when the lattice principal remaining is j-a in the up step
         //j - b in the down step 
         for(int a = 2; a >=0; a--){
            for(int b = 2; b >=0; b--){
            
            //moves onto next iteration of loop in exception cases

               //can't make paymemt of 10 less than 10 principal remaining 
               if (j < 2 && (b == 2 || a == 2))
                  continue; 

               //cannot make no payment if not ahead on payments 
               if ( j >= 40-n && (a == 0||b == 0))
                  continue; 
            
               //c1 is cash flow in the current period for the up step
               // c2 is cash flow in current period for the down step 
               //coupon * remaining principal + principal payment 
               c1 = coupon*5*j + 5*a;
               c2 = coupon*5*j + 5*b; 
            

               //calculates value of future cash flows 
               value = d[n][i] * ( 0.5 * (c1 + ValueFunkyBond(n+1, i+1,j-a)) + 
                              0.5* (c2 + ValueFunkyBond(n+1, i-1, j-b)));
            

               //gets the minimum value calculated thus far 
               if (value < min)
                  min = value; 
            

            }

         }
         //sets value of FunkyBond at (n,i,j) to be the minimum value found 
         //of all possible future cash flows 
         FunkyBond[n][i][j] = min; 
      }
   }

   

   //record that calculations have now been done at (n,i,j)
   I[n][i][j] = 1; 

   //returns value of funky bond
   return FunkyBond[n][i][j];

}

void AllocateMemory () {

   int n, i;

   // Allocate memory for term structure data.
   par = List(60);

   // Allocate memory for lattice arrays.
   d = AllocateLatticeArray ();
   V = AllocateLatticeArray ();

   // Allocate space for the Funky Bond valuation.
   // Indicator variable:
   I = (int ***) calloc (61, sizeof (int **));
   for (n = 0; n <= 60; n++) {
      I[n] = (int **) calloc (2*n+1, sizeof (int *));
      I[n] += n;
      for (i = -n; i <= n; i++) {
         I[n][i] = (int *) calloc (21, sizeof (int));
      }
   }
   // Valuation:
   FunkyBond = (double ***) calloc (61, sizeof (double **));
   for (n = 0; n <= 60; n++) {
      FunkyBond[n] = (double **) calloc (2*n+1, sizeof (double *));
      FunkyBond[n] += n;
      for (i = -n; i <= n; i++) {
         FunkyBond[n][i] = (double *) calloc (21, sizeof (double));
      }
   }

   return;

}


