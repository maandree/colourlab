/*
  Jargon:  luminosity, luminance, luma, brightness, value, lightness:
           The degree to which a colour can be described as pure white.
           
           saturation, colourfulness, excitation purity, purity, chroma, (chromaticness):
	   The degree to which a colour can be described as non-grey.
	   
	   chromacity:
	   The saturation and the hue of a colour, that is the colours described without
	   regard to its luminosity.
	   
	   hue:
	   The degree to which a colour can be described as similar to or different from
	   the primary colours, regardless of saturation or luminosity.
	   
	   opposite colour, complementary colour:
	   The colours which to add to another to get to 0 saturation, with same luminosity,
	   pure greys have them self as their complementary colour.
	   
	   elementary colours:
	   The colours used as reference to describe all colours.
	   
	   primary colours:
	   The elementary colours with full saturation.
	   
	   unique hues, principal hues:
	   The hue for the primary colours.
	   
	   gamut:
	   The set of colours that can be represented.
	   
	   
   Comparison:  Munsell:
                Adds purple as a prinicpal hue to make the hue wheel perceptually uniform

 */


/*
            sRGB to Y'UV
	    
	    ⎛Y'⎞   ⎛  0,299      0,587    0,114 ⎞ ⎛R / 255⎞
	    ⎜U ⎟ = ⎜−0,147313  −0,28886   0,436 ⎟ ⎜G / 255⎟
	    ⎝V ⎠   ⎝  0,615    −0,51499  −1.0001⎠ ⎝B / 255⎠

            Y'UV to sRGB
	    
	    ⎛R⎞       ⎛1      0      1,13973⎞ ⎛Y'⎞
	    ⎜G⎟ = 255 ⎜1  −0,39465  −0,58060⎟ ⎜U ⎟
	    ⎝B⎠       ⎝1   2,03211      0   ⎠ ⎝V ⎠
	    
	    sRGB to YDbDr
	    
	    ⎛Y ⎞   ⎛ 0,299   0,587  0,114⎞ ⎛R / 255⎞
	    ⎜Db⎟ = ⎜−0,450  −0,883  1,333⎟ ⎜G / 255⎟
	    ⎝Dr⎠   ⎝−1,333   1,116  0,217⎠ ⎝B / 255⎠
	    
	    YDbDr to sRGB
	    
	    ⎛R⎞       ⎛1   0,00002303716148   −0.525912630661865⎞ ⎛Y ⎞
	    ⎜G⎟ = 255 ⎜1  −0,129132898890509   0,267899328207599⎟ ⎜Db⎟
	    ⎝B⎠       ⎝1   0,64679059978955   −0,000079202543533⎠ ⎝Dr⎠
	    
	    sRGB to YIQ
	    
	    ⎛Y⎞   ⎛ 0,299      0,587      0,114  ⎞ ⎛R / 255⎞
	    ⎜I⎟ = ⎜0,959716  −0,274453  −0,321263⎟ ⎜G / 255⎟
	    ⎝Q⎠   ⎝0,211456  −0,522591  0,311135 ⎠ ⎝B / 255⎠
	    
	    YIQ to sRGB
	    
	    ⎛R⎞       ⎛1   0,9563   0,6210⎞ ⎛Y⎞
	    ⎜G⎟ = 255 ⎜1  −0,2721  −0,6474⎟ ⎜I⎟
	    ⎝B⎠       ⎝1  −1,1070   1,7046⎠ ⎝Q⎠
	    
	    Y'PbPr to Y'CbCr
	    
	    (Y', Cb, Cr) = (16, 128, 128) + (219, 224, 224) ∗ (Y', Pb, Pr)
	    
	    Y'CbCr to Y'PbPr
	    
	    (Y', Pb, Pr) = [(Y', Cb, Cr) − (16, 128, 128)] ÷ (219, 224, 224)
	    
	    sRGB to Y'PbPr
	    
	    ⎛Y'⎞   ⎛  0,299      0,587      0,114  ⎞ ⎛R / 255⎞
	    ⎜Pb⎟ = ⎜−0,168736  −0,331264     0,5   ⎟ ⎜G / 255⎟
	    ⎝Pr⎠   ⎝   0,5     −0,418688  −0,081312⎠ ⎝B / 255⎠
	    
	    sRGB to Y'PbPr
	    
	    ⎛R⎞       ⎛1  −0,0000012195   1,4019995887⎞ ⎛Y'⎞
	    ⎜G⎟ = 255 ⎜1  −0,3441356779  −0,7141361556⎟ ⎜Pb⎟
	    ⎝B⎠       ⎝1   1,7720000661   0,0000004063⎠ ⎝Pr⎠
 */


/**
 * @author  Mattias Andrée, <a href="mailto:maandree@kth.se">maandree@kth.se</a>
 */
public class Colour //TRY TO KEEP this class optimised for speed
{
    /**
     * Polynomial interpolation degree for hues
     */
    private static final int DEGREE = 4;
    
    
    /**
     * Empirical elementary colours: (yellow, yellow–red, red, red–blue, blue, blue–green, green, green–yellow) at 30 % white and 50 % chromacity.
     */
    private static final int[][] elementary = {{205, 101, 108}, {164, 110, 176}, { 36, 149, 190}, {  0, 169, 159},
					       { 50, 166, 121}, {156, 173,  81}, {204, 173,  71}, {218, 128,  77}};
    
    /**
     * The intensity of some empirical pure grey colours
     */
    private static final double[] greysW = { .10, .20, .30, .40, .50, .60, .70, .80, .90, .95, .97 };
    
    /**
     * The standard RGB red value of the empirical pure grey colours
     */
    private static final double[] greysR = { 25., 65.72, 94.38, 116.66, 136.1, 155.2, 175.68, 198.72, 225, 237, 243 };
    
    /**
     * The standard RGB green value of the empirical pure grey colours
     */
    private static final double[] greysG = greysR;
    
    /**
     * The standard RGB blue value of the empirical pure grey colours
     */
    private static final double[] greysB = greysR;
    
    
    
    //boundary check:  0 ≤ 2⋅|lum − 0.5| + sat ≤ 1 ∧ 0 ≤ sat, lum, white, black ≤ 1
    
    /**
     * Constructor
     * 
     * @param  lum  The luminosity of the colour [0, 1]
     * @param  sat  The saturation of the colour [0, 1]
     * @param  hue  The hue of the colour to interpolate: [0, 400[ gon
     */
    public Colour(final double lum, final double sat, final double hue)
    {
	this.saturation = sat;
	
	double h = hue;
	while (h >= 400.)  h -= 400.;
	while (h < 0.)     h += 400.;
	this.hue = h >= 400. ? 0. : h;
	
	final double[] m = ((lastLum == (this.luminosity = lum)) && (lastGrey != null)) ? lastGrey
	                   : (lastGrey = calculateLinearGrey(lastLum = lum));
	
	final int[][] fixed = new int[8][3];
	for (int i = 0; i < 8; i++)
	{
	    final double[] e = toLinear(elementary[i][0], elementary[i][1], elementary[i][2]);
	    
	    final double rk = 2. * (e[0] - m[0]);
	    final double gk = 2. * (e[1] - m[1]);
	    final double bk = 2. * (e[2] - m[2]);
	    
	    double r = m[0] + rk * sat;
	    double g = m[1] + gk * sat;
	    double b = m[2] + bk * sat;
	    
	    fixed[i] = toStandard(r, g, b);
	}
	
	final int[] srgb = hueToColour(hue, fixed);
	
	this.sRgb = srgb[0];
	this.srGb = srgb[1];
	this.srgB = srgb[2];
    }
    
    
    
    // Cache variables
    
    /**
     * Linear RGB red value coefficients for pure greys
     */
    private static double[] greyRk = null;
    
    /**
     * Linear RGB green value coefficients for pure greys
     */
    private static double[] greyGk = null;
    
    /**
     * Linear RGB blue value coefficients for pure greys
     */
    private static double[] greyBk = null;
    
    /**
     * The luminosity of the created colour
     */
    private static double lastLum = 0;
    
    /**
     * The zero chromacity corresponding colour, encoded in linear RGB, to the created colour
     */
    private static double[] lastGrey = null;
    
    
    // Display colour variables
    
    /**
     * Red value in the standard RGB model
     */
    private final int sRgb;
    
    /**
     * Green value in the standard RGB model
     */
    private final int srGb;
    
    /**
     * Blue value in the standard RGB model
     */
    private final int srgB;
    
    
    // Colour variables
    
    /**
     * The colour's luminosity [0, 1]
     */
    private final double luminosity;
    
    /**
     * The colour's saturation [0, 1]: the amount of colour, distance from grey
     */
    private final double saturation;
    
    /**
     * The colour's hue in [0, 400[ gon: 0 gon = 100 % yellow; 100 gon = 100 % red; 200 gon = 100 % blue; 300 gon = 100 % green
     */
    private final double hue;
    
    
    
    /**
     * Construct colour from the parameters: whiteness, saturation, hue
     * 
     * @param   white  The colour's whiteness [0, 1]
     * @param   sat    The colour's saturation [0, 1]: the amount of colour, distance from grey
     * @param   hue    The colour's hue in [0, 400[ gon: 0 gon = 100 % yellow; 100 gon = 100 % red; 200 gon = 100 % blue; 300 gon = 100 % green
     * @return         The colour with the given values
     */
    public Colour usingWSH(final double white, final double sat, final double hue)
    {
	final double black = 1. - white - sat;
	final double lum = (1. + white + black) / 2.;
	return new Colour(lum, sat, hue);
    }
    
    /**
     * Construct colour from the parameters: blackness, saturation, hue
     * 
     * @param   black  The colour's blackness [0, 1]
     * @param   sat    The colour's saturation [0, 1]: the amount of colour, distance from grey
     * @param   hue    The colour's hue in [0, 400[ gon: 0 gon = 100 % yellow; 100 gon = 100 % red; 200 gon = 100 % blue; 300 gon = 100 % green
     * @return         The colour with the given values
     */
    public Colour usingBSH(final double black, final double sat, final double hue)
    {
	final double white = 1. - black - sat;
	final double lum = (1. + white + black) / 2.;
	return new Colour(lum, sat, hue);
    }
    
    /**
     * Construct colour from the parameters: luminosity, saturation, hue
     * 
     * @param   lum  The colour's luminosity [0, 1]
     * @param   sat  The colour's saturation [0, 1]: the amount of colour, distance from grey
     * @param   hue  The colour's hue in [0, 400[ gon: 0 gon = 100 % yellow; 100 gon = 100 % red; 200 gon = 100 % blue; 300 gon = 100 % green
     * @return       The colour with the given values
     */
    public Colour usingLSH(final double lum, final double sat, final double hue)
    {
        //final double white = lum + sat / 2.;
	//final double black = 1. - white - sat;
	return new Colour(lum, sat, hue);
    }
    
    /**
     * Construct colour from the parameters: blackness, whiteness, hue
     * 
     * @param   black  The colour's blackness [0, 1]
     * @param   white  The colour's whiteness [0, 1]
     * @param   hue    The colour's hue in [0, 400[ gon: 0 gon = 100 % yellow; 100 gon = 100 % red; 200 gon = 100 % blue; 300 gon = 100 % green
     * @return         The colour with the given values
     */
    public Colour usingBWH(final double black, final double white, final double hue)
    {
	final double lum = (1. + white + black) / 2.;
	final double sat = black + white - 1;
	return new Colour(lum, sat, hue);
    }
    
    
    /**
     * Calculates a pure grey in linear RGB
     * 
     * @param   whiteness  The colour's whiteness or luminosity (or addative inverse blackness)
     * @return             The grey encoded in linear RGB
     */
    private static double[] calculateLinearGrey(final double whiteness)
    {
	final int n = greysW.length;
	if (greyRk == null)
	{
	    final double[][] x = new double[n][n];
	    for (int i = 0; i < n; i++)
	    {
		double e = 1, de = greysW[i];
		for (int j = 0; j < n; j++)
		{
		    x[i][j] = e;
		    e *= de;
		}
	    }
	    
	    final double[] gr = new double[n], gg = new double[n], gb = new double[n];
	    for (int i = 0; i < n; i++)
	    {
		final double[] rgb = toLinear(greysR[i], greysG[i], greysB[i]);
		gr[i] = rgb[0];
		gg[i] = rgb[1];
		gb[i] = rgb[2];
	    }
	    
	    greyRk = eliminate(x, gr);
	    greyGk = eliminate(x, gg);
	    greyBk = eliminate(x, gb);
	}
	
	double r = 0, g = 0, b = 0, e = 1;
	for (int i = 0; i < n; i++)
	{
	    r += greyRk[i] * e;
	    g += greyBk[i] * e;
	    b += greyGk[i] * e;
	    e *= whiteness;
	}
	
	return new double[] { r, g, b };
    }
    
    /**
     * Interpolates colour with a hue, using fixed colours evenly distributed in hue starting at 0
     * 
     * @param   hue    The hue of the colour to interpolate: [0, 400[ gon
     * @param   fixed  Fixed colours: {0 gon, 50 gon, 100 gon, ... 350 gon} × {red, green, blue}
     * @return         sRGB colour components: {red, green, blue}
     */
    private static int[] hueToColour(final double hue, final int[][] fixed)
    {
        double[][] f = new double[8][];
        double[] frk = new double[8];
        double[] fgk = new double[8];
	double[] fbk = new double[8];
	double[] fi;
	int[] fxi;
        for (int i = 0; i < 8; i++)
	{
	    frk[i] = (fi = f[i] = toLinear((fxi = fixed[i])[0], fxi[1], fxi[2]))[0];
	    fgk[i] = fi[1];
	    fbk[i] = fi[2];
	}
	
        int $hue = (int)(hue);
	$hue = ($hue % 400) + 400;
        if ($hue < 450)
            $hue += 400;
	int midH = ($hue / 50) & 7;
        int lowH = (midH - 1) & 7;
	
        int[] hs = new int[2];
        hs[0] = lowH & 7;
	hs[1] = (lowH + 1) & 7;
	
        double[][] xs = new double[DEGREE][DEGREE];
        double[] rk = new double[DEGREE];
	double[] gk = new double[DEGREE];
        double[] bk = new double[DEGREE];
        for (int y = 0; y < DEGREE; y++)
        {
	    int ym = (lowH + y) & 7;
	    rk[y] = frk[ym];
	    gk[y] = fgk[ym];
	    bk[y] = fbk[ym];
	    double c = 1, m = 400 + (lowH + y) * 50;
	    for (int x = 0; x < DEGREE; x++)
		xs[y][x] = c *= m;
	}
	rk = eliminate(xs, rk);
        gk = eliminate(xs, gk);
        bk = eliminate(xs, bk);
	
        double h = 1, r = 0, g = 0, b = 0;
	
        for (int i = 0; i < DEGREE; i++)
        {
	    r += (h *= $hue) * rk[i];
	    g +=  h          * gk[i];
	    b +=  h          * bk[i];
	}
	
        return toStandard(r, g, b);
    }
    
    
    /**
     * Gaussian elimination
     * 
     * @param   x  Square matrix
     * @param   y  Matrix augment
     * @return     Coefficients
     */
    private static double[] eliminate(final double[][] x, final double[] y)
    {
        int n = x.length;
        final double[] r = new double[n];
        final double[][] b = new double[n][n];
	
        System.arraycopy(y, 0, r, 0, n);
        for (int i = 0; i < n; i++)
            System.arraycopy(x[i], 0, b[i], 0, n);
	
        for (int k = 0, m = n - 1; k < m; k++)
            for (int i = k + 1; i < n; i++)
	    {
		double mul = b[i][k] / b[k][k];
		for (int j = k + 1; j < n; j++)
		    b[i][j] -= b[k][j] * mul;
		r[i] -= r[k] * mul;
	    }
	
        for (int k = n - 1; k > 0; k--)
            for (int i = 0; i < k; i++)
                r[i] -= r[k] * b[i][k] / b[k][k];
	
        for (int k = 0; k < n; k++)
            r[k] /= b[k][k];
	
        return r;
    }
    
    
    /**
     * Converts sRGB [0, 255] to linear RGB [0, 1]
     * 
     * @param   r  The red   intensity
     * @param   g  The green intensity
     * @param   b  The blue  intensity
     * @return     Linear RGB colours components
     */
    private static double[] toLinear(final int r, final int b, final int g)
    {
        return new double[] {
	            Math.pow(r / 255., 0.439764585),
		    Math.pow(g / 255., 0.439764585),
		    Math.pow(b / 255., 0.439764585)
	        };
    }
    
    
    /**
     * Converts sRGB [0, 255] to linear RGB [0, 1]
     * 
     * @param   r  The red   intensity
     * @param   g  The green intensity
     * @param   b  The blue  intensity
     * @return     Linear RGB colours components
     */
    private static double[] toLinear(final double r, final double b, final double g)
    {
        return new double[] {
	            Math.pow(r / 255., 0.439764585),
		    Math.pow(g / 255., 0.439764585),
		    Math.pow(b / 255., 0.439764585)
	        };
    }
    
    /**
     * Converts linear RGB [0, 1] to sRGB [0, 255]
     * 
     * @param   r  The red   intensity
     * @param   g  The green intensity
     * @param   b  The blue  intensity
     * @return     sRGB colours components
     */
    private static int[] toStandard(final double r, final double b, final double g)
    {
	return new int[] {
	            (int)(0.5 + 255. * Math.pow(r, 2.273943909)),
		    (int)(0.5 + 255. * Math.pow(g, 2.273943909)),
		    (int)(0.5 + 255. * Math.pow(b, 2.273943909))
		};
    }
    
    /**
     * Converts from sRGB to CIELAB
     * 
     * @param   red           The red   intensity [0, 255]
     * @param   green         The green intensity [0, 255]
     * @param   blue          The blue  intensity [0, 255]
     * @return                CIELAB colour components: {L*, a*, b*}
     */
    private static double[] toLab(final int red, final int green, final int blue)
    {
        int ir = red  ;  if (ir < 0)  ir += 1 << 8;
        int ig = green;  if (ig < 0)  ig += 1 << 8;
        int ib = blue ;  if (ib < 0)  ib += 1 << 8;
	
        double r = ir / 255.;  r = r <= 0.4045 ? r / 12.92 : Math.pow((r + 0.055) / 1.055, 2.4);
        double g = ig / 255.;  g = g <= 0.4045 ? g / 12.92 : Math.pow((g + 0.055) / 1.055, 2.4);
        double b = ib / 255.;  b = b <= 0.4045 ? b / 12.92 : Math.pow((b + 0.055) / 1.055, 2.4);
	
        double x = (0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / 0.95047;
        double y = (0.2126729 * r + 0.7151522 * g + 0.0721750 * b);
        double z = (0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / 1.08883;
	
        x = x > 0.00885642 ? Math.pow(x, 1. / 3.) : (7.78 + 703. / 99900.) * x + 0.1379310;
        y = y > 0.00885642 ? Math.pow(y, 1. / 3.) : (7.78 + 703. / 99900.) * y + 0.1379310;
        z = z > 0.00885642 ? Math.pow(z, 1. / 3.) : (7.78 + 703. / 99900.) * z + 0.1379310;
	
        double rcL = 116 * y - 16;
        double rca = 500 * (x - y);
        double rcb = 200 * (y - z);
	
        return new double[] {rcL, rca, rcb};
    }
    
    
    /**
     * Gets the sRGB representation of the colour<br/>
     * <b>Note the the colour may be out of gamut, in which case one of the values are outside [0, 255]</b>
     * 
     * @return  { red, green, blue }
     */
    public int[] getStandardRGB()
    {
	return new int[] { this.sRgb, this.srGb, this.srgB };
    }
    
    /**
     * Gets the linear RGB representation of the colour
     * 
     * @return  { red, green, blue }
     */
    public double[] getLinearRGB()
    {
	return toLinear(this.sRgb, this.srGb, this.srgB);
    }
    
    /**
     * Gets the CIELAB representation of the colour
     * 
     * @return  { L*, a*, b* }
     */
    public double[] getCIELAB()
    {
	return toLab(this.sRgb, this.srGb, this.srgB);
    }
    
    
    /**
     * Gets the colour's hun in gon: 0 gon = 100 % yellow; 100 gon = 100 % red; 200 gon = 100 % blue; 300 gon = 100 % green
     * 
     * @return  The colour's hue in [0, 400[ gon
     */
    public double getHue()
    {
	return this.hue;
    }
    
    /**
     * Gets the colour's saturation: the amount of colour, distance from grey
     * 
     * @return  The colour's saturation [0, 1]
     */
    public double getSaturation()
    {
	return this.saturation;
    }
    
    /**
     * Gets the colour's luminosity: the brightness, where only 0.5 can reach full chromacity
     * 
     * @return  The colour's luminosity [0, 1]
     */
    public double getLuminosity()
    {
	return this.luminosity;
    }
    
    /**
     * Gets the colour's whiteness
     * 
     * @return  The colour's whiteness [0, 1]
     */
    public double getWhiteness()
    {
	return this.luminosity + this.saturation / 2.;
    }
    
    /**
     * Gets the colour's blackness
     * 
     * @return  The colour's blackness [0, 1]
     */
    public double getBlackness()
    {
	return 1. - this.getWhiteness() - this.saturation;
    }
    
    /**
     * Gets the colour's yellowness
     * 
     * @return  The colour's yellowness [0, 1]
     */
    public double getYellowness()
    {
	if ((  0. <= this.hue) && (this.hue < 100.))  return (100. - this.hue) / 100.;
	if ((300. <= this.hue) && (this.hue < 400.))  return (this.hue - 300.) / 100.;
	return 0.;
    }
    
    /**
     * Gets the colour's redness
     * 
     * @return  The colour's redness [0, 1]
     */
    public double getRedness()
    {
	if ((  0. <= this.hue) && (this.hue < 100.))  return this.hue / 100.;
	if ((100. <= this.hue) && (this.hue < 200.))  return (200. - this.hue) / 100.;
	return 0.;
    }
    
    /**
     * Gets the colour's blueness
     * 
     * @return  The colour's blueness [0, 1]
     */
    public double getBlueness()
    {
	if ((100. <= this.hue) && (this.hue < 200.))  return (this.hue - 100.) / 100.;
	if ((200. <= this.hue) && (this.hue < 300.))  return (300. - this.hue) / 100.;
	return 0.;
    }
    
    /**
     * Gets the colour's greenness
     * 
     * @return  The colour's greenness [0, 1]
     */
    public double getGreenness()
    {
	if ((200. <= this.hue) && (this.hue < 300.))  return (this.hue - 200.) / 100.;
	if ((300. <= this.hue) && (this.hue < 400.))  return (400. - this.hue) / 100.;
	return 0.;
    }
    
}

