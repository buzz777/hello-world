#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14159265359;

double f(double, double);

double Euler(double, double, double);

double Predictor(double, double, double, double);

double Runge_Kutt_3(double, double, double);

double Runge_Kutt_4(double, double, double);

double ABM(double, double, double, double, double, double);

double Gir(double, double, double, double, double, double);

double eps(double, double);

double abs_d(double);

int main()
{
    int N = 400;
    double a = 0;
    double b = pi;
    double h = (b - a) /(double) N;
    double y, y1, y2, y3, y_new = 0;
    double eps_curr, eps_total = 0;
    double lambda = 0.000000001;
    int i = 0;

    FILE* f_out = fopen("output_Euler.txt", "w");

/**EULER**/

    y = cos(a);
    fprintf(f_out, "%d\n", N);
    for (i = 0; i < N; i++){
        y = Euler(a+h*i, y, h);
        eps_curr = eps(a+h*i, y);

        fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a + h*i), eps_curr);
        if (eps_curr > eps_total){
            eps_total = eps_curr;
        }
    }
    fprintf(f_out, "\n\neps total: %.15lf\nh: %.15lf\n", eps_total, h);
    fclose(f_out);

/**PREDICTOR**/
    eps_total = 0;
    y = cos(a);
    f_out = fopen("output_Predictor.txt", "w");
    fprintf(f_out, "%d\n", N);
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a), eps_curr);
    for (i = 1; i < N; i++){
        y = Predictor(a+h*(i-1), y, h, lambda);
        eps_curr = eps(a+h*i, y);
        fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a + h*i), eps_curr);
        if (eps_curr > eps_total)
            eps_total = eps_curr;
    }
    fprintf(f_out, "\n\neps total: %.15lf\nh: %.15lf\n", eps_total, h);
    fclose(f_out);

/**RUNGE-KUTT 3RD**/
    eps_total = 0;
    y = cos(a);
    f_out = fopen("output_Runge-Kutt_3.txt", "w");
    fprintf(f_out, "%d\n", N);
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a), eps_curr);
    for (i = 1; i < N; i++){
        y = Runge_Kutt_3(a+h*(i-1), y, h);
        eps_curr = eps(a+h*i, y);
        fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a + h*i), eps_curr);
        if (eps_curr > eps_total)
            eps_total = eps_curr;
    }
    fprintf(f_out, "\n\neps total: %.15f\nh: %.15lf\n", eps_total, h);
    fclose(f_out);

/**RUNGE-KUTT 4TH**/
    eps_total = 0;
    y = cos(a);
    f_out = fopen("output_Runge-Kutt_4.txt", "w");
    fprintf(f_out, "%d\n", N);
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a), eps_curr);
    for (i = 1; i < N; i++){
        y = Runge_Kutt_4(a+h*(i-1), y, h);
        eps_curr = eps(a+h*i, y);
        fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a + h*i), eps_curr);
        if (eps_curr > eps_total)
            eps_total = eps_curr;
    }
    fprintf(f_out, "\n\neps total: %.15f\nh: %.15lf\n", eps_total, h);
    fclose(f_out);

/**ABM**/
    eps_total = 0;

    y3 = cos(a);
    y2 = Runge_Kutt_4(a, y3, h);
    y1 = Runge_Kutt_4(a + h, y2, h);
    y = Runge_Kutt_4(a + 2*h, y1, h);

    f_out = fopen("output_ABM.txt", "w");
    fprintf(f_out, "%d\n", N);
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y3, cos(a), eps(a, y3));
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y2, cos(a+h), eps(a+h, y2));
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y1, cos(a+2*h), eps(a+2*h, y1));
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a+3*h), eps(a+3*h, y));
    for (i = 4; i < N; i++){
        y_new = ABM(a+h*(i-1), y, y1, y2, y3, h);
        y3 = y2;
        y2 = y1;
        y1 = y;
        y = y_new;
        eps_curr = eps(a+h*i, y);
        fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a + h*i), eps_curr);
        if (eps_curr > eps_total)
            eps_total = eps_curr;
    }
    fprintf(f_out, "\n\neps total: %.15f\nh: %.15lf\n", eps_total, h);
    fclose(f_out);

/**GIR**/
    eps_total = 0;

    y3 = cos(a);
    y2 = Runge_Kutt_4(a, y3, h);
    y1 = Runge_Kutt_4(a + h, y2, h);
    y = Runge_Kutt_4(a + 2*h, y1, h);


    f_out = fopen("output_GIR.txt", "w");
    fprintf(f_out, "%d\n", N);
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y3, cos(a), eps(a, y3));
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y2, cos(a+h), eps(a+h, y2));
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y1, cos(a+2*h), eps(a+2*h, y1));
    fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a+3*h), eps(a+3*h, y));
    for (i = 4; i < N; i++){
        y_new = Gir(a+h*(i), y, y1, y2, y3, h);
        y3 = y2;
        y2 = y1;
        y1 = y;
        y = y_new;
        eps_curr = eps(a+h*i, y);
        fprintf(f_out, "%.15lf %.15lf %.15lf\n", y, cos(a + h*i), eps_curr);
        if (eps_curr > eps_total)
            eps_total = eps_curr;
    }
    fprintf(f_out, "\n\neps total: %.15f\nh: %.15lf\n", eps_total, h);
    fclose(f_out);

    return 0;
}

double f(double x, double y){
    double k = 2.5;
    return -sin(x)+k*(y - cos(x));
}

double Euler(double x, double y, double h){
    return (y + h*f(x, y));
}

double Predictor(double x, double y, double h, double lambda){
    double y_curr_s = 0;
    double y_curr_s1 = 2;
    int s = 0;

    while (s < 10000 && (abs_d(y_curr_s - y_curr_s1)/abs_d(y_curr_s1)) > lambda){
        if(s == 0){
            y_curr_s = Euler(x, y, h);
        }
        else{
            y_curr_s1 = y_curr_s;
            y_curr_s = y + h*((f(x, y) + f(x+h, y_curr_s1))/2.0);
        }
        s++;
    }

    return y_curr_s;
}

double Runge_Kutt_3(double x, double y, double h){
    double result = 0;
    double k1, k2, k3 = 0;

    k1 = h*f(x, y);
    k2 = h*f(x + h/2.0,y + k1/2.0);
    k3 = h*f(x + h, y - k1 + 2*k2);
    result = y + (k1 + 4*k2 + k3)/6.0;

    return result;
}

double Runge_Kutt_4(double x, double y, double h){
    double result = 0;
    double k1, k2, k3, k4 = 0;

    k1 = h*f(x, y);
    k2 = h*f(x + h/2.0,y + k1/2.0);
    k3 = h*f(x + h/2.0,y + k2/2.0);
    k4 = h*f(x + h, y + k3);
    result = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;

    return result;
}

double ABM(double x, double y, double y1, double y2, double y3, double h){
    double result = 0;
    double y_w = y + h*(55*f(x,y) - 59*f(x - h, y1) + 37*f(x - 2*h, y2) - 9*f(x - 3*h, y3))/24.0;

    result = y + h*(9*f(x+h, y_w) + 19*f(x, y) - 5*f(x - h, y1) + f(x - 2*h, y2))/24.0;

    return result;
}

double Gir(double x, double y, double y1, double y2, double y3, double h){
    double result = 0;
    double k = 2.5;

    result = (48*y - 36*y1 + 16*y2 - 3*y3 - 12*h*(sin(x) + k*cos(x)))/(double)(25 - 12*h*k);

    return result;
}

double eps(double x, double y){
    return abs_d(cos(x) - y);
}

double abs_d(double x){
    if(x >= 0)
        return x;
    return -x;
}
