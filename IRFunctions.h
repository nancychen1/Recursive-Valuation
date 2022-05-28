

void     Calibrate ();
double **AllocateLatticeArray ();
double   Bond (int, double);
double   ParRate (int, int, int);
void     MakeParCurves ();
double   Hump (double);
double   Twist (double);
void     OutputFiles (void);
void     ComputeCurves ();

////////////////////////////////////////////////////////////////////////////////
// This function calibrates the state-dependent single-period discount
//    factors d[n][i] to explain the current par yield curve in par[].
////////////////////////////////////////////////////////////////////////////////
void Calibrate () {

   int n, i;
   double lower, upper, trial, gamma;

   // Specify the gamma parameter with delta_t = 0.5.
   gamma = exp (sigma * sqrt(0.5));

   // Calibrate going forward in time.
   for (n = 0; n < 60; n++) {

      // Bounds on rho.
      lower = 0.0;
      upper = 1.0;

      // Iterate using interval bisection until close.
      while (1) {

         // Compute trial discount factors at period n for this trial rho.
         trial = (lower + upper) / 2.0;
         for (i = -n; i <= n; i += 2) {
            d[n][i] = pow (trial, pow(gamma, i));
         }

         // Value a bond maturing in period n+1 using the trial rho.
         Bond (n+1, par[n+1]);

         // See if error tolerance has been met. If so, terminate search.
         if (fabs (V[0][0] - 100.0) < 0.00001) {
            break;
         }

         else if (V[0][0] > 100.0) {
            // Too little discounting -- trial is too big. The true rho is
            //    less than trial so trial becomes an upper bound for rho.
            upper = trial;
         }

         else {
            // Too much discounting.
            lower = trial;
         }

      }

   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function values a non-callable bond that matures at period m and pays
//    an annual coupon of "coupon" using the current d[n][i]'s for 0 <= n < m.
////////////////////////////////////////////////////////////////////////////////
double Bond (int m, double coupon) {

   int n, i;
   double C;

   // Boundary conditions.
   for (i = -m; i <= m; i += 2) {
      V[m][i] = 0.0;
   }

   // Make coupon semiannual.
   coupon /= 2.0;

   // Iterate backward recursively.
   for (n = m-1; n >= 0; n--) {
      for (i = -n; i <= n; i += 2) {

         // Calculate the cash flow.
         if (n+1 == m) {
            C = coupon + 100.0;   // Return principal at maturity.
         } else {
            C = coupon;
         }

         V[n][i] = d[n][i] * (0.5 * (C + V[n+1][i-1]) + 0.5 * (C + V[n+1][i+1]));

      }
   }

   return V[0][0];

}

////////////////////////////////////////////////////////////////////////////////
// This function calculates the par rate at lattice node (n,i)
//    expressed as a percent for a bond that matures at period m.
// Must have 0 <= n < m <= 60, and -n <= i <= n. n and i must
//    have the same parity (odd/even).
////////////////////////////////////////////////////////////////////////////////
double ParRate (int n, int i, int m) {

   double C, v0, v1;

   // Value a zero coupon bond maturing at time m.
   Bond (m, 0.0);
   v0 = V[n][i];

   // Value a bond paying a coupon of $1 per year.
   Bond (m, 1.0);
   v1 = V[n][i];


   // The fair coupon C at node (n,i) for a par bond maturing in period m
   //    satisfies: C * v1 + (1-C) * v0 = 100.0; so...

   // General case:
   if (v1-v0 > 1e-100) {
      C = (100.0 - v0) / (v1 - v0);
   }

   // To avoid overflow (this happens near the top of the lattice where rates are huge):
   else {
      C = 1e+100;   // This is so rare that any number here will work.
   }

   return C;

}




////////////////////////////////////////////////////////////////////////////////
// Allocate memory for a lattice array x[n][i], where 0 <= n <= 60 and
//    -n <= i <= n.
////////////////////////////////////////////////////////////////////////////////
double **AllocateLatticeArray () {

   int n;
   double **x;

   x = (double **) calloc (61, sizeof (double *));
   for (n = 0; n <= 60; n++) {
      x[n] = (double *) calloc (2*n + 1, sizeof (double));
      x[n] += n;
   }

   return x;

}

////////////////////////////////////////////////////////////////////////////////
// This function creates a par yield curve in percent.
// The curve is reported to the array "par[]".
////////////////////////////////////////////////////////////////////////////////
void MakeParCurve () {

   int n;
   double A, B, C, t, shortrate, longrate, hump;
   FILE *fp;

   // "shortrate" is the rate at the short end of the yield curve (time 0).
   // "longrate"  is the rate at the long  end of the yield curve (time 30).
   // "hump", drags the curve up by that amount at time 4 years,
   //   leaving the short and long rates fixed.

   printf ("Please specify the par yield curve in percent...\n");
   shortrate = GetDouble ("What is the short rate?... ");
   longrate  = GetDouble ("What is the long rate?.... ");
   hump      = GetDouble ("What is the hump?......... ");

   // Solve for the parameters A and B to satisfy:
   //   A + B * Twist (0.5)  + C * Hump (0.5)  = "shortrate",  and
   //   A + B * Twist (30.0) + C * Hump (30.0) = "longrate", when C = "hump".
   C = hump;
   B = (shortrate - longrate + C * (Hump (30.0) - Hump (0.5)))
       / (Twist (0.5) - Twist (30.0));
   A = shortrate - B * Twist (0.5) - C * Hump (0.5);

   // Step forward in time period by period.
   for (n = 1; n <= 60; n++) {

      // Maturity in years rather than periods.
      t = 0.5 * n;

      // Par rate to time t, as determined by the parameters.
      par[n] = A + B * Twist(t) + C * Hump(t);

   }

   // We now have par[1] = shortrate, and par[60] = longrate.

   // Report "shortrate" and "longrate" for viewing.
   fp = fopen ("DataPoints.txt", "w");
   fprintf (fp, "%6.2f %8.5f  %6.2f %8.5f  %6.2f %8.5f\n",
                0.5, shortrate, 4.0, par[8], 30.0, longrate);
   fclose (fp);


   return;

}

////////////////////////////////////////////////////////////////////////////////
// Nelson-Seigel Hump function.
////////////////////////////////////////////////////////////////////////////////
double Hump (double t) {

   double y, tau = 4.0;

   // The hump function ((1.0 - exp(-t)) / t) - exp(-t) maximizes at t = 1.79328
   //    and its maximum value is 0.298426.

   // Scale time so that the hump maximizes at t = tau years.
   t *= 1.79328 / tau;

   // Scale "y" so that the maximum value is 1.0    (1.0 / 0.298426 = 3.35092).
   if (t == 0) {
      y = 0.0;
   } else {
      y = 3.35092 * ( ((1.0 - exp(-t)) / t) - exp(-t) );
   }

   return (y);

}

////////////////////////////////////////////////////////////////////////////////
// Twist function.
////////////////////////////////////////////////////////////////////////////////
double Twist (double t) {

   double y, tau = 4.0;

   // Scale time as in the "Hump" function.
   t *= 1.79328 / tau;

   if (t == 0) {
      y = 1.0;
   } else {
      y = (1.0 - exp(-t)) / t;
   }

   return (y);

}

////////////////////////////////////////////////////////////////////////////////
// This function generates some .txt files for viewing data with Curves.tex.
// The par curve must be constructed and the lattice must be calibrated.
////////////////////////////////////////////////////////////////////////////////
void OutputFiles () {

   int i, n;
   double t, r, v0, v1;
   FILE *fp;

   // Report the lattice rate data for viewing with TeX software.
   // Show only for short rates r <= 15%, with semiannual compounding.
   // This corresponds to d[n][i] >= (1.0 + 0.15/2)^(-1) = 0.9302326.
   fp = fopen ("LatticeRates.txt", "w");
   for (n = 0; n < 60; n++) {
      t= 0.5 * n;
      for (i = -n; i <= n; i += 2) {
         if (d[n][i] > 0.9302326) {
            r = ParRate (n, i, n+1);
            fprintf (fp, "%8.2f %8.2f\n", t, r);
         }
      }
   }
   fclose (fp);

   // Report the par curve.
   fp = fopen ("par.txt", "w");
   for (n = 1; n <= 60; n++) {
      t = 0.5 * n;
      fprintf (fp, "%4.1f %8.6f\n", t, par[n]);
   }
   fclose (fp);

   // Report zero coupon curve.
   fp = fopen ("zero.txt", "w");
   for (n = 1; n <= 60; n++) {
      t = 0.5 * n;
      v0 = Bond (n, 0.0) / 100.0;
      r = 200.0 * (pow(v0, -1.0/n) - 1.0);
      fprintf (fp, "%4.1f %8.6f\n", t, r);
   }
   fclose (fp);

   // Report single-period forward curve.
   fp = fopen ("forward.txt", "w");
   for (n = 0; n < 60; n++) {
      t = 0.5 * n;
      if (n == 0) v0 = 100.0;
      else        v0 = Bond (n, 0.0);
      v1 = Bond (n+1, 0.0);
      r =  200.0 * ((v0/v1 - 1.0));
      fprintf (fp, "%4.1f %8.6f\n", t, r);
   }
   fclose (fp);

   // Report "shortrate", "longrate", and the 4-year rate for viewing.
   fp = fopen ("DataPoints.txt", "w");
   fprintf (fp, "%6.2f %8.5f  %6.2f %8.5f  %6.2f %8.5f\n",
                0.5, par[1], 4.0, par[8], 30.0, par[60]);
   fclose (fp);


   return;

}

////////////////////////////////////////////////////////////////////////////////
// From the par curve, this function computes and reports the spot curve and
//    the single-period forward rates.
////////////////////////////////////////////////////////////////////////////////
void ComputeCurves () {

   int i, n;
   double coupons, *discount, zero_n, forward_n;
   FILE *fp;

   // Allocate memory for the discount function.
   discount = List (60);

   // Report the par curve.
   fp = fopen ("par.txt", "w");
   for (n = 1; n <= 60; n++) {
      fprintf (fp, "%4.1f %8.6f\n", 0.5*n, par[n]);
   }
   fclose (fp);

   // Step 1: compute the discount function from the par curve.

   discount[0] = 1.0;
   for (n = 1; n <= 60; n++) {
      coupons = 0.0;
      for (i = 1; i < n; i++) {
         coupons += par[n] * 0.5 * discount[i];
      }
      discount[n] = (100.0 - coupons) / (100.0 + par[n] * 0.5);
   }

   // Step 2: now that the discount function discount[] is computed...

   // ... compute and report the zero coupon curve,
   fp = fopen ("zero.txt", "w");
   for (n = 1; n <= 60; n++) {
      zero_n = 200.0 * (pow(discount[n], -1.0 / n) - 1.0);
      fprintf (fp, "%4.1f %8.6f\n", 0.5*n, zero_n);
   }
   fclose (fp);

   // ... and finally compute and report single-period forward rates.
   fp = fopen ("forward.txt", "w");
   for (n = 0; n < 60; n++) {
      forward_n = 200.0 * (pow (discount[n+1]/discount[n], -1) - 1.0);
      fprintf (fp, "%4.1f %8.6f\n", 0.5*n, forward_n);
   }
   fclose (fp);

   free (discount);

}




