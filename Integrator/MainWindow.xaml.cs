using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Integrator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void TheButton_Click(object sender, RoutedEventArgs e)
        {
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            double right = 0;
            double left = 0;
            double trap = 0;
           // x.Add(0);
           // x.Add(0);
            //x.Add(1);
            x.Add(1);
            x.Add(2);
            //x.Add(3);
            // y.Add(2);
            //y.Add(-2);
            //y.Add(2);
            y.Add(-1);
            y.Add(1);
           // y.Add(2);
            Integration.Curves.MonotonicallyIncreasing c = new Integration.Curves.MonotonicallyIncreasing(x, y);
            if (c.Validate().Errors.Count() == 0)
            {
                left = c.RiemannSum(false, .00001);
                System.Diagnostics.Debug.Print(left.ToString());// lefthand
                //System.Diagnostics.Debug.Print(S)
                for(int i = 0; i<= 1000000; i++)
                {
                    right = c.RiemannSum(true, .00001);
                }
                
                System.Diagnostics.Debug.Print(right.ToString());//righthand 
                //System.Diagnostics.Debug.Print(c.RiemannSum(true,1.25,1.75).ToString());// righthand
                //System.Diagnostics.Debug.Print(c.RiemannSum(false,1.25,1.75).ToString());// lefthand
                trap = c.TrapazoidalRiemannSum();
                System.Diagnostics.Debug.Print(trap.ToString());// Trapazoidal
                double rightprecision = (right - trap) / trap;
                System.Diagnostics.Debug.Print(rightprecision.ToString());
                double leftprecision = (left - trap) / trap;
                System.Diagnostics.Debug.Print(leftprecision.ToString());
                //System.Diagnostics.Debug.Print(c.TrapazoidalRiemannSum(1.25,1.75).ToString());// Trapazoidal
            }
            else
            {
                System.Diagnostics.Debug.Print("bad news bears!");
            }

        }
    }
}
