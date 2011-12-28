import javax.swing.*;
import java.awt.event.*;

public class UserInterface extends JFrame implements ActionListener
{
	static UserInterface ui;
	
	public static void main( String args[] )
	{
		try
		{
//			UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
		} catch (Exception e)
		{
			e.printStackTrace(System.out);
		}
		ui = new UserInterface();
	}

	public void actionPerformed(ActionEvent e) {
		
	}
	
}
