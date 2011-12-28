
public class Example extends  LinearSlider
{

	public static void main(String[] args) 
	{
		double steps, pos;

		steps = getSteps();
		System.out.println("Current step: " + steps);
		
		pos = getPos();
		System.out.println("Current dist: " + pos );
		
		setHome(1000);
		
		steps = getSteps();
		System.out.println("Current step: " + steps);
		
		pos = getPos();
		System.out.println("Current dist: " + pos );
		
		
		
	}
	
	public static void CallingStaticClass()
	{

		steps = LinearSlider.getSteps();
		System.out.println("Current step: " + steps);
		
		pos = LinearSlider.getPos();
		System.out.println("Current dist: " + pos );
		
		LinearSlider.setHome(1000);
		
		steps = LinearSlider.getSteps();
		System.out.println("Current step: " + steps);
		
		pos = LinearSlider.getPos();
		System.out.println("Current dist: " + pos );
		
		
		
	}

}
