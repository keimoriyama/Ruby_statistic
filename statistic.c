#include <ruby/ruby.h>
#include <stdio.h>
#include <stdlib.h>

void check_type(VALUE elem){
    switch(TYPE(elem)){
        case T_FLOAT:
            break;
        case T_BIGNUM:
            break;
        case T_FIXNUM:
            break;
        default:
            rb_raise(rb_eTypeError, "not valid value");
            break;
    }
}

static VALUE
ary_mean(VALUE self)
{
    double e, sum = 0;
    VALUE ary = rb_ary_to_ary(self);
    int len = RARRAY_LENINT(ary);
    int i;
    VALUE elem;
    for(i = 0; i < len; i++){
       elem = RARRAY_AREF(ary, i);
       check_type(elem);
       e = NUM2DBL(elem);
       sum = sum + e;
    }
    sum = sum / len;
    return DBL2NUM(sum);
}

VALUE static
rb_ary_var(VALUE self){
    VALUE e;
    double elem, variance, t_val, cmean, l;
    int i;
    variance = 0;
    VALUE mean = ary_mean(self);
    cmean = NUM2DBL(mean);
    VALUE ary = self;
    l = RARRAY_LENINT(ary);
    for (i = 0; i < l; i++) {
        e = RARRAY_AREF(ary, i);   
        check_type(e);
        elem = NUM2DBL(e);
        t_val = elem-cmean;
        variance = variance + t_val*t_val;
        /* printf("%f\n", variance); */
    }
    variance = variance / l;
    return DBL2NUM(variance);
}

#include <math.h>

double fact(int n){
    double value = 1;
    while(n >= 0){
        value = value * n;
        n--;
    }
    return value;
}

double beta(double p, double q){
    double value = 0;
    value = tgamma(p)*tgamma(q)/tgamma(p+q);
    return value;
}

//t分布の密度関数
double t_densty_func(double t, int df){
    return pow((1+t*t/df), -1.0*(df+1)/2)/(sqrt(df)*beta(df/2, 0.5));
}

double t_prob(double end_p, int df){
    double start_t = 0;
    double diff = 0.01;
    double p=0;
    end_p /= 2;
    while(p < end_p) {
        p = p + (t_densty_func(start_t, df)+t_densty_func(start_t+diff, df))*diff/2;
        start_t = start_t + diff;
    }
    p = 2*p;
    return start_t;
}

double kai_double(double x){
    if(x == 0)return 0;
    else
    return exp(-x/2)/(sqrt(2*x)*tgamma(0.5));
}

double kai_double_p(double p){
    double start = 0;
    double value = 0;
    double step = 0.0001;
    while(value < p){
        value = value + ((kai_double(start+step) + kai_double(start))*step)/2;
        start = start + step;
    }
    return start;
}

VALUE static
rb_kai_double(VALUE self, VALUE rb_ary2, VALUE alpha){
    VALUE rb_ary1 = rb_ary_to_ary(self);
    double a = NUM2DBL(alpha);
    if(TYPE(rb_ary2) != T_ARRAY || a > 1){
        rb_raise(rb_eTypeError, "not valid value");
    }
    double table[2][2], f[2][2];
    table[0][0] = NUM2DBL(RARRAY_AREF(rb_ary1, 0));
    table[0][1] = NUM2DBL(RARRAY_AREF(rb_ary1, 1));
    table[1][0] = NUM2DBL(RARRAY_AREF(rb_ary2, 0));
    table[1][1] = NUM2DBL(RARRAY_AREF(rb_ary2, 1));
    /* return DBL2NUM(kai_double_p(0, 10)); */
    double sum = 0;
    for(int i = 0; i < 2; i++){
        for(int j =0; j < 2; j++){
            sum = sum + table[i][j];
        }
    }
    f[0][0] = (table[0][0] + table[0][1])*(table[0][0]+table[1][0])/sum;
    f[0][1] = (table[0][0] + table[0][1])*(table[0][1]+table[1][1])/sum;
    f[1][0] = (table[0][0] + table[1][0])*(table[1][0] + table[1][1])/sum;
    f[1][1] = (table[0][1] + table[1][1])*(table[1][1] + table[1][0])/sum;
    double kai = 0, tmp;
    for(int i = 0; i < 2; i++){
        for(int j =0; j < 2; j++){
           tmp = table[i][j] - f[i][j];
           kai = kai + tmp*tmp/f[i][j];
        }
    }
    double low_x = kai_double_p(1-a);
    int res = 0;
    if(low_x <= kai)
       res = 1;
    else
        res = -1;
    printf("カイ二乗値:%.3f\n", kai);
    printf("有意水準:%.2f\n", a);
    printf("棄却域：kai>%.3f\n", low_x);
    return INT2NUM(res);
}

VALUE static
rb_ary_t_mesure(VALUE self, VALUE rb_mu, VALUE rb_alpha){
    VALUE rb_mean = ary_mean(self);
    VALUE rb_var = rb_ary_var(self);
    double alpha = NUM2DBL(rb_alpha);
    if(alpha > 1 || alpha < 0)
        rb_raise(rb_eTypeError, "not valid value");
    double l = RARRAY_LEN(self);
    double df = l - 1.0;
    double mu = NUM2DBL(rb_mu);
    double mean = NUM2DBL(rb_mean);
    double var = NUM2DBL(rb_var);
    double t_fr = mean - mu;
    double t_ds = sqrt(var/l);
    // t統計量
    double t_value = t_fr/t_ds;
    // 棄却域
    double p = t_prob(1-alpha, df);
    printf("t統計量：%.3f\n", t_value);
    printf("有意水準:%.2f\n", alpha);
    printf("棄却域：|t|>%.3f\n", p);
    int res=0;
    t_value = fabs(t_value);
    if(p < t_value)
        res = 1;
    else
        res = -1;
    return INT2NUM(res);
}

VALUE static
rb_ary_std(VALUE self){
    VALUE rb_var = rb_ary_var(self);
    double var = NUM2DBL(rb_var);
    double std = sqrt(var);
    return DBL2NUM(std);
}

void
Init_statistic(void){
    rb_define_method(rb_cArray, "mean",ary_mean, 0);
    rb_define_method(rb_cArray, "var",rb_ary_var, 0);
    rb_define_method(rb_cArray, "std",rb_ary_std, 0);
    rb_define_method(rb_cArray, "t_mesure",rb_ary_t_mesure, 2);
    rb_define_method(rb_cArray, "kai_mesure",rb_kai_double, 2);
}
