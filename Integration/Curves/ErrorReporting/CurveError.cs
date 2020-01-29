using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Integration.Curves.ErrorReporting
{
    public class CurveError
    {
        #region Properties
        public string Message { get; set; }
        public int Row { get; set; }
        public int Column { get; set; }
        #endregion

        #region Constructors
        #endregion
        public CurveError(string message, int row, int col)
        {
            Message = message;
            Row = row;
            Column = col;
        }
    }
}
