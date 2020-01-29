using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Integration.Curves
{
    public abstract class PairedData
    {
        #region Notes
        #endregion
        #region Fields

        #endregion
        #region Properties
        public double[] X;
        public double[] Y;
        public bool IsFunctional{ get; set;}
        #endregion
        #region Constructors

        #endregion
        #region Voids

        #endregion
        #region Functions
        public abstract double GetXFromY(double y);
        public abstract double GetYFromX(double x);
        public abstract ErrorReporting.ErrorReport Validate();
        public double RiemannSum(bool UseRightHandSums)
        {

            //https://en.wikipedia.org/wiki/Riemann_sum
            int n = X.Count();
            double result = 0;
            if (!UseRightHandSums)
            {
                for (int i = 0; i <= n - 2; i++)
                {
                    result += (Y[i]) * (X[i+1] - X[i]);
                }
            }
            else
            {
                for (int i = 1; i <= n-1; i++)
                {
                    result += (Y[i]) * (X[i] - X[i-1]);
                }
            }                
            return result;
        }
        public double RiemannSum(bool UseRightHandSums, double precisionTarget)
        {

            //https://en.wikipedia.org/wiki/Riemann_sum
            int n = X.Count();
            //precision must be pro rated by ordinate. makes sense to pro rate it by slope.
            double perOrdinatePrecision = precisionTarget;// precisionTarget / (X.Count() - 1);
            double result = 0;
            double trapazoidalArea = 0;
            double absolutePrecision = 0;
            if (!UseRightHandSums)
            {
                for (int i = 0; i <= n - 2; i++)
                {
                    double dx = (X[i + 1] - X[i]);
                    double y = Y[i];
                    //calculate actual area
                    trapazoidalArea = (Y[i + 1] + y) / 2 * (dx);
                    //calculate num slices based on slope
                    double slope = Slope(i, i + 1);
                    
                    int slices = 1;
                    if (slope != 0)//if slope equals 0, the error reserved for that segment could be reassigned to the remainder of the function.
                    {
                        absolutePrecision = trapazoidalArea * perOrdinatePrecision;
                        slices = N(slope, dx, absolutePrecision);
                    }
                    
                    // calculate the y values for each dx
                    
                    double ddx = dx / slices;
                    
                    for (int j = 0; j < slices; j++)
                    {
                        
                        result += y * ddx;
                        y += (slope * ddx);
                    }
                }
            }
            else
            {
                for (int i = 1; i <= n - 1; i++)
                {
                    //calculate actual area
                    double dx = (X[i] - X[i-1]);
                    double y = Y[i - 1];// + (slope * ddx);
                    trapazoidalArea = (Y[i] + y) / 2 * (dx);

                    //calculate num slices based on slope
                    double slope = Slope(i-1, i);
                    
                    int slices = 1;
                    if(slope != 0)
                    {
                        absolutePrecision = trapazoidalArea * perOrdinatePrecision;
                        slices = N(slope, dx, absolutePrecision);
                    }
                    
                    // calculate the y values for each dx
                    
                    double ddx = dx / slices;
                    
                    for (int j = 0; j < slices; j++)
                    {
                        y += (slope * ddx);
                        result += y * ddx;
                    }
                }
            }
            return result;
        }
        private double Slope(int idx1, int idx2)
        {
            double top = Y[idx1] - Y[idx2];
            double bottom = X[idx1] - X[idx2];
            return  top/ bottom;
        }
        private int N(double slope, double dx, double precision)
        {
            double rsqrd = dx * dx;
            double n = Math.Ceiling((slope * rsqrd) / (2 * precision));
            double absn = Math.Abs(n);
            return (int)absn;

        }
        public double RiemannSum(bool UseRightHandSums, double precision, double xStart, double xEnd)
        {
            if (xStart == xEnd) { return 0; }
            if (xEnd < X[0]) { throw new ArithmeticException("The range of the paired data relationship is greater than the range requested"); }
            if (xStart > X.Last()) { throw new ArithmeticException("The range of the paired data relationship is less than the range requested"); }
            if (xStart > xEnd) { throw new ArithmeticException("The range start is greater than the range end"); }// -result; } //requires looping through the ordinates backwards
            if (xStart < X[0]) { throw new ArithmeticException("A portion of the range is outside of the range of the function"); }
            if (xEnd > X.Last()) { throw new ArithmeticException("A portion of the range is outside of the range of the function"); }
            double result = 0;
            double tmpX = 0;
            double slope = 0;
            double trapazoidalArea = 0;
            int n = X.Count();
            bool hasValues = false;
            int largestXindexNotInInterval = -1;
            int smallestXindexNotInInterval = -1;
            double y = 0;
            double absolutePrecision = 0;
            double perOrdinatePrecision = precision;
            if (!UseRightHandSums)
            {
                for (int i = 0; i <= n - 2; i++)
                {
                    tmpX = X[i];
                    if (tmpX >= xStart)
                    {
                        if (tmpX <= xEnd)
                        {
                            double dx = (X[i + 1] - X[i]);
                            y = Y[i];
                            //calculate actual area
                            trapazoidalArea = (Y[i + 1] + y) / 2 * (dx);
                            //calculate num slices based on slope
                            slope = Slope(i, i + 1);

                            int slices = 1;
                            if (slope != 0)//if slope equals 0, the error reserved for that segment could be reassigned to the remainder of the function.
                            {
                                absolutePrecision = trapazoidalArea * perOrdinatePrecision;
                                slices = N(slope, dx, absolutePrecision);
                            }

                            // calculate the y values for each dx

                            double ddx = dx / slices;

                            for (int j = 0; j < slices; j++)
                            {

                                result += y * ddx;
                                y += (slope * ddx);
                            }
                            hasValues = true;
                        }
                        else
                        {
                            largestXindexNotInInterval = i;
                            break;
                        }
                    }
                    else
                    {
                        smallestXindexNotInInterval = i;
                    }
                }

                if (!hasValues)
                {
                    //inbetween ordinates
                    //evaluate from left side (xStart) of bin?
                    if (smallestXindexNotInInterval == -1)
                    {
                        smallestXindexNotInInterval = 0;
                    }
                    if (largestXindexNotInInterval == -1)
                    {
                        largestXindexNotInInterval = smallestXindexNotInInterval + 1;
                    }
                    double dx = xEnd - xStart;
                    slope = ((Y[largestXindexNotInInterval] - Y[smallestXindexNotInInterval]) / (X[largestXindexNotInInterval] - X[smallestXindexNotInInterval]));
                    y = Y[smallestXindexNotInInterval] + (slope * (xStart - X[smallestXindexNotInInterval]));
                    return y * dx;
                }
                else
                {
                    //add the first rectangle and the last rectangle.
                    slope = 0;
                    if (largestXindexNotInInterval == -1)
                    {
                        //then the range is between the last and second to last ordinates.
                        if (xEnd == X.Last())
                        {
                            //dont add anything.
                        }
                        else
                        {
                            result += Y[n - 2] * (xEnd - X[n - 2]);
                        }
                    }
                    else
                    {
                        result += Y[largestXindexNotInInterval - 1] * (xEnd - X[largestXindexNotInInterval - 1]);
                    }
                    if (smallestXindexNotInInterval == -1)
                    {
                        //then the range is between the first two ordinates
                        if (xStart == X[0])
                        {
                            //dont add anything
                            return result;
                        }
                        else
                        {
                            slope = (Y[1] - Y[0]) / (X[1] - X[0]);
                            y = (Y[0] + (slope * (xStart - X[0])));//need to figure out the Y value for the xStart
                            result += y * (X[1] - xStart);
                        }
                    }
                    else
                    {
                        slope = (Y[smallestXindexNotInInterval + 1] - Y[smallestXindexNotInInterval]) / (X[smallestXindexNotInInterval + 1] - X[smallestXindexNotInInterval]);
                        y = (Y[smallestXindexNotInInterval] + ((slope) * (xStart - X[smallestXindexNotInInterval])));//need to figure out the Y value fo rthe xStart
                        result += y * (X[smallestXindexNotInInterval + 1] - xStart);
                    }
                }
            }
            else
            {
                for (int i = 1; i <= n - 1; i++)
                {
                    tmpX = X[i];
                    if (tmpX >= xStart)
                    {
                        if (tmpX <= xEnd)
                        {
                            //calculate actual area
                            double dx = (X[i] - X[i - 1]);
                            y = Y[i - 1];// + (slope * ddx);
                            trapazoidalArea = (Y[i] + y) / 2 * (dx);

                            //calculate num slices based on slope
                            slope = Slope(i - 1, i);

                            int slices = 1;
                            if (slope != 0)
                            {
                                absolutePrecision = trapazoidalArea * perOrdinatePrecision;
                                slices = N(slope, dx, absolutePrecision);
                            }

                            // calculate the y values for each dx

                            double ddx = dx / slices;

                            for (int j = 0; j < slices; j++)
                            {
                                y += (slope * ddx);
                                result += y * ddx;
                            }
                            hasValues = true;
                        }
                        else
                        {
                            largestXindexNotInInterval = i;
                            break;
                        }
                    }
                    else
                    {
                        smallestXindexNotInInterval = i;
                    }
                }

                if (!hasValues)
                {
                    //inbetween ordinates
                    //evaluate from right side (xEnd) of bin?
                    if (smallestXindexNotInInterval == -1)
                    {
                        smallestXindexNotInInterval = 0;
                    }
                    if (largestXindexNotInInterval == -1)
                    {
                        largestXindexNotInInterval = smallestXindexNotInInterval + 1;
                    }
                    double dx = xEnd - xStart;
                    y = Y[smallestXindexNotInInterval] + ((Y[largestXindexNotInInterval] - Y[smallestXindexNotInInterval]) / (X[largestXindexNotInInterval] - X[smallestXindexNotInInterval]) * (xEnd - X[smallestXindexNotInInterval]));
                    return y * dx;
                }
                else
                {
                    //add the first rectangle and the last rectangle.

                    if (largestXindexNotInInterval == -1)
                    {
                        //then the range is between the last and second to last ordinates.
                        if (xEnd == X.Last())
                        {
                            //dont add anything.
                        }
                        else
                        {
                            // calculate the value of y at xEnd
                            slope = (Y[n - 1] - Y[n - 2]) / (X[n - 1] - X[n - 2]);
                            y = (Y[n - 2] + (slope * (xEnd - X[n - 2])));//need to figure out the Y value for the xEnd
                            result += Y[n - 2] * (xEnd - X[n - 2]);
                        }
                    }
                    else
                    {
                        result += Y[largestXindexNotInInterval - 1] * (xEnd - X[largestXindexNotInInterval - 1]);
                    }
                    if (smallestXindexNotInInterval == -1)
                    {
                        //then the range is between the first two ordinates
                        if (xStart == X[0])
                        {
                            //dont add anything
                            return result;
                        }
                        else
                        {

                            result += Y[1] * (X[1] - xStart);
                        }
                    }
                    else
                    {
                        slope = (Y[smallestXindexNotInInterval + 1] - Y[smallestXindexNotInInterval]) / (X[smallestXindexNotInInterval + 1] - X[smallestXindexNotInInterval]);
                        y = (Y[smallestXindexNotInInterval] + ((slope) * (X[smallestXindexNotInInterval]) - xStart));//need to figure out the Y value fo rthe xStart
                        result += Y[smallestXindexNotInInterval + 1] * (X[smallestXindexNotInInterval + 1] - xStart);
                    }
                }
            }

            return result;
        }
        public double RiemannSum(bool UseRightHandSums, double xStart, double xEnd)
        {
            if (xStart == xEnd) { return 0; }
            if (xEnd < X[0]) { throw new ArithmeticException("The range of the paired data relationship is greater than the range requested"); }
            if (xStart > X.Last()) { throw new ArithmeticException("The range of the paired data relationship is less than the range requested"); }
            if (xStart > xEnd) { throw new ArithmeticException("The range start is greater than the range end"); }// -result; } //requires looping through the ordinates backwards
            if (xStart < X[0]) { throw new ArithmeticException("A portion of the range is outside of the range of the function"); }
            if (xEnd > X.Last()) { throw new ArithmeticException("A portion of the range is outside of the range of the function"); }
            double result = 0;
            double tmpX = 0;
            double slope = 0;
            int n = X.Count();
            bool hasValues = false;
            int largestXindexNotInInterval = -1;
            int smallestXindexNotInInterval = -1;
            double y = 0;
            if (!UseRightHandSums)
            {
                for (int i = 0; i <=n-2 ; i++)
                {
                    tmpX = X[i];
                    if (tmpX >= xStart)
                    {
                        if (tmpX <= xEnd)
                        {
                            result += (Y[i])*(X[i + 1] - X[i]);
                            hasValues = true;
                        }
                        else
                        {
                            largestXindexNotInInterval = i;
                            break;
                        }
                    }
                    else
                    {
                        smallestXindexNotInInterval = i;
                    }
                }

                    if (!hasValues)
                    {
                    //inbetween ordinates
                    //evaluate from left side (xStart) of bin?
                    if (smallestXindexNotInInterval == -1)
                    {
                        smallestXindexNotInInterval = 0;
                    }
                    if(largestXindexNotInInterval == -1)
                    {
                        largestXindexNotInInterval = smallestXindexNotInInterval+1;
                    }
                    double dx = xEnd - xStart;
                    slope = ((Y[largestXindexNotInInterval] - Y[smallestXindexNotInInterval]) / (X[largestXindexNotInInterval] - X[smallestXindexNotInInterval]));
                        y = Y[smallestXindexNotInInterval] +  (slope* (xStart-X[smallestXindexNotInInterval]));
                        return y * dx;
                }
                else
                {
                    //add the first rectangle and the last rectangle.
                    slope = 0;
                    if (largestXindexNotInInterval == -1)
                    {
                        //then the range is between the last and second to last ordinates.
                        if (xEnd == X.Last())
                        {
                            //dont add anything.
                        }
                        else
                        {
                            result += Y[n-2]* (xEnd - X[n-2]);
                        }
                    }
                    else
                    {
                        result += Y[largestXindexNotInInterval-1]* (xEnd - X[largestXindexNotInInterval-1]);
                    }
                    if (smallestXindexNotInInterval == -1)
                    {
                        //then the range is between the first two ordinates
                        if (xStart == X[0])
                        {
                            //dont add anything
                            return result;
                        }
                        else
                        {
                            slope = (Y[1] - Y[0]) / (X[1] - X[0]);
                            y = (Y[0] + (slope * (xStart - X[0])));//need to figure out the Y value for the xStart
                            result += y * (X[1] - xStart);
                        }
                    }
                    else
                    {
                        slope = (Y[smallestXindexNotInInterval + 1] - Y[smallestXindexNotInInterval]) / (X[smallestXindexNotInInterval + 1] - X[smallestXindexNotInInterval]);
                        y = (Y[smallestXindexNotInInterval] + ((slope) * (xStart-X[smallestXindexNotInInterval])));//need to figure out the Y value fo rthe xStart
                        result += y * (X[smallestXindexNotInInterval+1]-xStart);
                    }
                }
            }else
            {
                for (int i = 1; i <= n - 1; i++)
                {
                    tmpX = X[i];
                    if (tmpX >= xStart)
                    {
                        if (tmpX <= xEnd)
                        {
                            result += (Y[i]) * (X[i] - X[i-1]);
                            hasValues = true;
                        }
                        else
                        {
                            largestXindexNotInInterval = i;
                            break;
                        }
                    }
                    else
                    {
                        smallestXindexNotInInterval = i;
                    }
                }

                if (!hasValues)
                {
                    //inbetween ordinates
                    //evaluate from right side (xEnd) of bin?
                    if (smallestXindexNotInInterval == -1)
                    {
                        smallestXindexNotInInterval = 0;
                    }
                    if (largestXindexNotInInterval == -1)
                    {
                        largestXindexNotInInterval = smallestXindexNotInInterval + 1;
                    }
                    double dx = xEnd - xStart;
                    y = Y[smallestXindexNotInInterval] + ((Y[largestXindexNotInInterval] - Y[smallestXindexNotInInterval]) / (X[largestXindexNotInInterval] - X[smallestXindexNotInInterval]) * (xEnd- X[smallestXindexNotInInterval]));
                    return y * dx;
                }
                else
                {
                    //add the first rectangle and the last rectangle.
                    
                    if (largestXindexNotInInterval == -1)
                    {
                        //then the range is between the last and second to last ordinates.
                        if (xEnd == X.Last())
                        {
                            //dont add anything.
                        }
                        else
                        {
                            // calculate the value of y at xEnd
                            slope = (Y[n-1] - Y[n-2]) / (X[n-1] - X[n-2]);
                            y = (Y[n-2] + (slope * (xEnd - X[n-2])));//need to figure out the Y value for the xEnd
                            result += Y[n - 2] * (xEnd - X[n - 2]);
                        }
                    }
                    else
                    {
                        result += Y[largestXindexNotInInterval - 1] * (xEnd - X[largestXindexNotInInterval - 1]);
                    }
                    if (smallestXindexNotInInterval == -1)
                    {
                        //then the range is between the first two ordinates
                        if (xStart == X[0])
                        {
                            //dont add anything
                            return result;
                        }
                        else
                        {

                            result += Y[1] * (X[1] - xStart);
                        }
                    }
                    else
                    {
                        slope = (Y[smallestXindexNotInInterval + 1] - Y[smallestXindexNotInInterval]) / (X[smallestXindexNotInInterval + 1] - X[smallestXindexNotInInterval]);
                        y = (Y[smallestXindexNotInInterval] + ((slope) * (X[smallestXindexNotInInterval]) - xStart));//need to figure out the Y value fo rthe xStart
                        result += Y[smallestXindexNotInInterval+1] * (X[smallestXindexNotInInterval+1] - xStart);
                    }
                }
            }

            return result;
        }
        public double TrapazoidalRiemannSum()
        {
            int iter = X.Count() - 2;
            double result = 0;
            for(int i = 0; i <= iter; i++)
            {
                result += (Y[i + 1] + Y[i]) / 2 * (X[i + 1] - X[i]);
            }
            return result;
        }
        public double TrapazoidalRiemannSum(double xStart, double xEnd)
        {
            if (xStart == xEnd) { return 0; }
            if (xEnd < X[0]) { throw new ArithmeticException("The range of the paired data relationship is greater than the range requested"); }
            if (xStart > X.Last()) { throw new ArithmeticException("The range of the paired data relationship is less than the range requested"); }
            if (xStart > xEnd){ throw new ArithmeticException("The range start is greater than the range end"); }// -result; } //requires looping through the ordinates backwards
            if (xStart < X[0]) { throw new ArithmeticException("A portion of the range is outside of the range of the function"); }
            if (xEnd > X.Last()) { throw new ArithmeticException("A portion of the range is outside of the range of the function"); }
            int iter = X.Count() - 2;
            double result = 0;
            double tmpX = 0;
            bool hasValues = false;
            int largestXindexNotInInterval = -1;
            int smallestXindexNotInInterval = -1;
            for (int i = 0; i <= iter; i++)
            {
                tmpX = X[i];
                if (tmpX >= xStart)
                {
                    if(tmpX<= xEnd)
                    {
                        result += (Y[i + 1] + Y[i]) / 2 * (X[i + 1] - X[i]);
                        hasValues = true;
                    }else
                    {
                        largestXindexNotInInterval = i;
                        break;
                    }
                }
                else
                {
                    smallestXindexNotInInterval = i;
                }
            }

            if (!hasValues)
            {
                //inbetween ordinates
                //slope must be constant
                if (smallestXindexNotInInterval == -1)
                {
                    smallestXindexNotInInterval = 0;
                }
                if (largestXindexNotInInterval == -1)
                {
                    largestXindexNotInInterval = smallestXindexNotInInterval + 1;
                }
                double dx = xEnd - xStart;
                double top = Y[largestXindexNotInInterval] + Y[smallestXindexNotInInterval];
                return top / 2 * dx;
            }else
            {
                //add the first trapazoid and the last trapazoid.
                //only evaluate on original range of the function.
                double frontPartTop = 0;
                double slope = 0;
                if (largestXindexNotInInterval == -1) {
                        //then the range is between the last and second to last ordinates.
                    if (xEnd == X.Last())
                    {
                        //dont add anything.
                    }else
                    {
                        slope = (Y[iter + 1] - Y[iter]) / (X[iter + 1] - X[iter]);
                        frontPartTop = Y[iter+1] + (Y[iter] + ((slope) * (xEnd - X[iter])));//need to figure out the y value for the xEnd...
                        result += frontPartTop / 2 * (xEnd - X[iter]);
                    }
                }else
                {
                    slope = (Y[largestXindexNotInInterval] - Y[largestXindexNotInInterval-1]) / (X[largestXindexNotInInterval] - X[largestXindexNotInInterval-1]);
                    frontPartTop = Y[largestXindexNotInInterval] + (Y[largestXindexNotInInterval-1]+((slope)*(xEnd-X[largestXindexNotInInterval])));//need to figure out the y value for the xEnd...
                    result += frontPartTop / 2 * (xEnd - X[largestXindexNotInInterval-1]);
                }
                if (smallestXindexNotInInterval == -1) {
                    //then the range is between the first two ordinates
                    if (xStart == X[0])
                    {
                        //dont add anything
                        return result;
                    }
                    else
                    {
                        slope = (Y[1] - Y[0]) / (X[1] - X[0]);
                        frontPartTop = Y[1] + (Y[0]+(slope*(xStart - X[0])));//need to figure out the Y value fo rthe xStart
                        result += frontPartTop / 2 * (X[1] - xStart);
                    }
                }
                else
                {
                    slope = (Y[smallestXindexNotInInterval + 1] - Y[smallestXindexNotInInterval]) / (X[smallestXindexNotInInterval + 1] - X[smallestXindexNotInInterval]);
                    frontPartTop = Y[smallestXindexNotInInterval+1] + (Y[smallestXindexNotInInterval]+((slope)* (xEnd - X[smallestXindexNotInInterval])));//need to figure out the Y value fo rthe xStart
                    result += frontPartTop / 2 * (X[1] - xStart);
                }
            }
            
            return result;
        }
        #endregion
    }
}
