# The 1 comp. PK model with dose-dependent bioavailability.

# [INPUT]
ka_pop = logit(0.2);
V_pop = log(20);
Cl_pop = 0.2;

beta_V_AGE = 1;
beta_ka_1 = 1;
AGE_med = 50;

omega_ka = 1;
omega_V = 1;
omega_Cl = 1;

Cc_a = 1;
Cc_b = 1;

ka_tv = expit(ka_pop);
V_tv = exp(V_pop);
Cl_tv = Cl_pop;

ka = ka_tv * expit(omega_ka) * (beta_ka_1 * SEX);
V = V_tv * exp(omega_V) * (AGE/AGE_med)^beta_V_AGE;
Cl = Cl_tv + omega_Cl

# [MODEL]
# PK model definition

### Priors
MW = 627.57;  					# g/mol for INCB161734

### Explicit functions
Cc = Ac/V;
Cc_nM = Cc/MW*1000000;			# nmol/L

### Initial conditions
Ad(0) = 0;
Ac(0) = 0;

### Differential equations
d/dt(Ad) = - ka*Ad;
d/dt(Ac) = ka*Ad - Cl*Ac/V;

Cc_Res_err = Cc * (1 + Cc_b ) + Cc_a;
