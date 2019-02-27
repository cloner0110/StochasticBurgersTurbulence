function x1 = SBE_FBM(beta,m)
 
x           = sqrt(m)*randn(1,m);
L           = m/2 ;
freqM       = [ 1 1:L (L-1):-1:1 ];
xft         = fft(x);
xft(1)      = 0; xft(L+1) = 0;  %Set the Nyquist to zero
x1ft        = xft .* (freqM .^ ( beta/2 ) );
x1          = real( ifft(x1ft) );
