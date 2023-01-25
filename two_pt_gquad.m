function approx = two_pt_gquad(f, a, b)
    approx = f( (b+a)/2 - (1/sqrt(3))*((b-a)/2) ) + f( (b+a)/2 + (1/sqrt(3))*((b-a)/2) ) ;