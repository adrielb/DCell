
public class LinearSlider
{
	// conversion factor in steps per inch 
	final static double FACTOR = 80000./5.;

	private static SerialControl sc = 
		new SerialControl("LinearSlider", "COM6",'/');
	
	// height z is in inches
	public static void move( double z )
	{
		move( z * FACTOR );
	}
	
	// height z in microsteps
	public static void move( int z )
	{
		sc.write("/1A" + z + "R\r");
	}
	
	// move in relative up direction
	public static void moveUp( double z )
	{
		moveUp( z * FACTOR );
	}
	
	public static void moveUp( int z )
	{
		if( z == 0 ) return;
		sc.write("/1P" + z + "R\r");
	}
	
	// move down relative to current position
	public static void moveDown( double z )
	{
		moveDown( z * FACTOR );
	}
	
	public static void moveDown( int z )
	{
		if( z == 0 ) return;
		sc.write("/1D" + z + "R\r");
	}
	
	// set zero location
	public static void setHome( double z )
	{
		setHome( z * FACTOR );
	}
	
	public static void setHome( int z )
	{
		sc.write("/1z"+ z + "R\r");
	}
	
	// Terminate current command or loop
	public static void terminate()
	{
		sc.write( "/1T" );
	}

	// Get position in steps
	public static int getSteps()
	{
		sc.write("/1?0");
		Integer.parseInt(sc.getAvailableData().trim());
		return pos;
	}
	
	// get position in inches
	public static double getPos()
	{
		return getSteps() / FACTOR;
	}
	
	
}
