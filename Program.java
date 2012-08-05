import javax.swing.*;
import java.awt.*;
import java.awt.image.*;


/**
 * @author  Mattias Andrée, <a href="mailto:maandree@kth.se">maandree@kth.se</a>
 */
public class Program extends JFrame
{
    public Program()
    {
	super("Free perceptual colour model");
	this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	this.pack();
	final Insets in = this.getInsets();
	this.setSize(new Dimension(451 + in.left + in.right, 451 + in.top + in.bottom));
    }
    
    
    
    public static void main(final String... args)
    {
	for (int il = 0, ln = imgs.length - 1; il <= ln; il++)
	{
	    System.err.println("At image " + il);
	    final BufferedImage img = imgs[il] = new BufferedImage(401, 401, BufferedImage.TYPE_INT_ARGB);
	    final Graphics2D gg = img.createGraphics();
	    
	    gg.setColor(Color.WHITE);
	    gg.fillRect(0, 0, 401, 401);
	    
	    for (int y = -200; y <= 200; y++)
		for (int x = -200; x <= 200; x++)
		{
		    final double fy = y / 200.;
		    final double fx = x / 200.;
		    
		    if (fy * fy + fx * fx > 1.)
			continue;
		    
		    final double l = il / (double)ln,
			         s = Math.sqrt(fy * fy + fx * fx),
			         h = ((Math.atan2(y, x) * 200. / Math.PI) + 100.) % 400.;
		    
		    final double ł = (l - .5 < 0.) ? (.5 - l) : (l - .5);
		    if (2. * ł + s > 1.)
			continue;
		    
		    final Colour colour = new Colour(l, s, Colour.toNaturalHue(h));
		    final int[] rgb = colour.getStandardRGB();
		    int r = rgb[0],
			g = rgb[1],
			b = rgb[2];
		    
		    if ((r > 255) || (g > 255) || (b > 255) || (r < 0) || (g < 0) || (b < 0))
		    {
			r = r > 255 ? 255 : r < 0 ? 0 : r;
			g = g > 255 ? 255 : g < 0 ? 0 : g;
			b = b > 255 ? 255 : b < 0 ? 0 : b;
			r = 255 - r;
			g = 255 - g;
			b = 255 - b;
		    }
		    
		    gg.setColor(new Color(r, g, b));
		    gg.fillRect(200 + x, 200 + y, 1, 1);
		}
	}
	
	final Program frame = new Program();
	frame.setVisible(true);
	
	final Thread thread = new Thread()
	    {
		@Override
		public void run()
		{
		    int diff = 1;
		    try
		    {
			for (;;)
			{
			    Thread.sleep(5000 / (Program.imgs.length - 1));
			    Program.lum += diff;
			    if (Program.lum >= Program.imgs.length - 1)
			    {
				Program.lum = Program.imgs.length - 1;
			        diff = -diff;
			    }
			    else if (Program.lum <= 0)
			    {
				Program.lum = 0;
			        diff = -diff;
			    }
			    frame.repaint();
			}
		    }
		    catch (final Throwable err)
		    {}
		}
	    };
	thread.setDaemon(true);
	thread.start();
    }
    
    
    
    private static BufferedImage[] imgs = new BufferedImage[1 + 4];
    
    static int lum = 0;
    
    
    
    @Override
    public void paint(final Graphics gp)
    {
	final Insets in = this.getInsets();
	gp.drawImage(imgs[lum], 25 + in.left, 25 + in.top, null);
    }
    
}

