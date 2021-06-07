import("stdfaust.lib");


//---------------------------`bubble`--------------------------
// bubble(f0, trig) : produces a water drop bubble sound
//
// #### Usage
//
// ```
// bubble(f0, trig) : _
// ```
//
// Where:
//
// * ` f0 `: base frequency of bubble sound
// * `trig`: trigs the bubble sound on the rising front
//
// #### Example
//
// ```
// button("drop") : bubble(600) : _
// ```
//
// #### Reference:
//
// <http://www.cs.ubc.ca/~kvdoel/publications/tap05.pdf>
//------------------------------------------------------------
//Shifter taken from ollie's wicked OWL patches
// Copywrite "Oli Larkin (contact@olilarkin.co.uk)";
// IIR hilbert transform Emmanuel Favreau (via Miller Puckette)
// fastest
hilbertef(x) = real(x), imag(x)
with {
  biquad(a1,a2,b0,b1,b2) = + ~ conv2(a1,a2) : conv3(b0,b1,b2) 
  with {
    conv3(k0,k1,k2,x) = k0*x + k1*x' + k2*x'';
    conv2(k0,k1,x) = k0*x + k1*x';
  };
  real = biquad(-0.02569, 0.260502, -0.260502, 0.02569, 1) 
        : biquad(1.8685, -0.870686, 0.870686, -1.8685, 1) ;
  imag = biquad(1.94632, -0.94657, 0.94657, -1.94632, 1) 
      : biquad(0.83774, -0.06338, 0.06338,  -0.83774, 1) ;
};

freqshift(x, shift) = negative(x), positive(x)
with {
  negative(x) = real(x)*cosv - imag(x)*sinv;
  positive(x) = real(x)*cosv + imag(x)*sinv;
  real(x) = hilbert(x) : _ , !;
  imag(x) = hilbert(x) : ! , _;
  
  phasor(x) = fmod((x/float(ma.SR) : (+ : ma.decimal) ~ _), 1.)  * (ma.PI * 2);

  sinv = sin(phasor(shift));
  cosv = cos(phasor(shift));

  hilbert = hilbertef;
};

ssb(shift, x) = freqshift(x, shift) : _ , !; // only take one sideband



shift_amount = shift*shift_scalar;

string(f,d) = (+~(de.fdelay4(960,ma.SR/f-1) : _ <: _,_ ' :> / (2) : * (d) ) ) ;
gain = hslider("gain",0.22,0,1,0.01);
base = hslider("Res Base",50,10,120,1):si.smoo;
freq1 = hslider("Res freq1",0,0,24,1):si.smoo;
freq2 = hslider("Res freq2",3,0,24,1):si.smoo;
freq3 = hslider("Res freq3",7,0,24,1):si.smoo;
freq4 = hslider("Res freq4",9,0,24,1):si.smoo;
damp = hslider("Res damp",0.98,0,0.99,0.01):si.smoo;

invol = hslider("invol",0.98,0,0.99,0.01) :si.smoo;
cutoffHigh = hslider("cutoffHi",20,20,20000,1):si.smoo;
cutoffLow = hslider("cutoffLow",20000,20,20000,1):si.smoo;

l1 = hslider("Res l1",0.22,0,0.22,0.01);
l2 = hslider("Res l2",0.22,0,0.22,0.01);
l3 = hslider("Res l3",0.22,0,0.22,0.01);

msec = ma.SR/1000.0;
shiftl = hslider("Pitch Shift L ", 12, -24, +24, 0.1);
shiftr = hslider("Pitch Shift R ", 12, -24, +24, 0.1);
ws = hslider("Pitch Shift Window Size ", 50, 20, 1000, 1) * msec : si.smooth(ba.tau2pole(0.005));
mixP = hslider("PitchShiftMix", 0.5, 0, 1, 0.01) : si.smooth(ba.tau2pole(0.005));

xf = 20 * msec;

transpose (w, x, s, sig) = de.fdelay(65536, d,sig)*ma.fmin(d/x,1) + de.fdelay(65536,d+w,sig)*(1-ma.fmin(d/x,1))
with {
  i = 1 - pow(2, s/12);
  d = i : (+ : +(w) : fmod(_,w)) ~ _;
};
	
pitchShifter(l,r) = l,r <: *(1-mixP), *(1-mixP), transpose(ws, xf, shiftl, l)*mixP, transpose(ws, xf, shiftr, r)*mixP :> _,_;

mixD = hslider("DelayMix [knob:4]", 1.0, 0, 1, 0.01) : si.smooth(ba.tau2pole(0.005));
echof(l) =  l<: *(1-mixD),  +~(de.fdelay4(maxDelLength,dLength-1) : *(dfeedback))*mixD :> _;

duration = hslider("Duration [knob:1]", 200,1,200,1)*0.001 :  si.smoo;
dfeedback = hslider("Feedback [knob:2]", 0.78, 0 , 0.999 ,0.01);//
//damping =  hslider("Damping", 0.99,0,1,0.01);
filter = _ <: _,_ ' :> / (2);
maxDelLength = 12000;
dLength = ma.SR*duration;

nbits = hslider("BitCrush BITS",24,1,24,1);
scaler = float(2^nbits-1);
round(x) = floor(x+0.5);
bitcrusher(nbits,x) = x :abs: *(scaler) : round : /(scaler) * (2*(x>0)-1.0);

redux = hslider("BitCrush SamplFrq",192000,44100/64,192000,1); 
sampleRedux = ba.downSample(redux); 

comb = hgroup("Comb[0]", +~(@(delLength-1) : *(feedback)));

delLength = hslider("Duration[1]", 1,1,200,1);
feedback = hslider("Feedback[2]", 0.0, 0 , 0.85 ,0.01);


  

bubble(f0,trig) = os.osc(f) * (exp(-damp*time) : si.smooth(0.99))
	with {
		damp = 0.043*f0 + 0.0014*f0^(3/2);
		f = f0*(1+sigma*time);
		sigma = eta * damp;
		eta = 0.075;
		time = 0 : (select2(trig>trig'):+(1)) ~ _ : ba.samp2sec;
	};

vol = hslider("volume [unit:dB]", 0, -96, 0, 0.1) : ba.db2linear : si.smoo;
freq = hslider("freq [unit:Hz]", 1000, 20, 24000, 1);
shift = hslider("shift [unit:hz]", 2.2, -20., 20, 0.001);
shift_scalar = hslider("shift_scalar", 6., 1., 100, 0.1);
lr_offset = hslider("lr_offset", 0., 0., 1., 0.00001);
mix = hslider("mix",0.5,0,1,0.01) : si.smooth(ba.tau2pole(0.005));

process = vgroup("Oscillator", os.osc(500.0) * vol <: *(1-mix) , *(1-mix), ssb(shift_amount ,_ )*mix , ssb(shift_amount +lr_offset,_ )*mix :> (_,_) );

shifter(l, r) = l, r <: *(1-mix) , *(1-mix), ssb(shift_amount ,l )*mix , ssb(shift_amount +lr_offset,r )*mix :> _,_ ;

// process = ba.pulsen(1, 10000) : pm.djembe(60, 0.3, 0.4, 1) <: string(ba.midikey2hz(base+freq1),damp)*l1 , string(ba.midikey2hz(base+freq2),damp)*l2 , string(ba.midikey2hz(base+freq3),damp)*l3 :> _ <: shifter(_,_) : echof,echof;
// shifter(l, r) = l, r <: *(1-mix) , *(1-mix), ssb(shift_amount ,l )*mix , ssb(shift_amount +lr_offset,r )*mix :> _,_ ; 