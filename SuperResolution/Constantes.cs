using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SuperResolution
{
    public class Constantes
    {
        public static int MINARG = 5;
        public static double BREAK_NUMBER = 0.0;
        public static int DEFAULT_SF = 22050;               /* Hz. Sampling Frequency */
        public static double DEFAULT_SHIFT = 5.0;           /* ms */
        public static double DEFAULT_MIN_PITCH = 40.0;      /* Hz */
        public static double DEFAULT_MAX_PITCH = 400.0;     /* Hz */
        public static int DEFAULT_DECIMATION = 4;           /* samples */
        public static int DEFAULT_TSILENT = 120;            /* max. abs sample amplitude of noise */
        public static double DEFAULT_THIGH = 0.88;
        public static double DEFAULT_TMIN = 0.75;
        public static double DEFAULT_TMAX_RATIO = 0.85;
        public static double DEFAULT_TDH = 0.77;
        public static int UNVOICED = 0;                     /* segment classifications */
        public static int VOICED = 1;
        public static int SILENT = 2;
        public static int HOLD = 1;
        public static int HELD = 1;
        public static int SEND = 2;
        public static double SENT = 2;

        public static int BEGINNING = 1;
        public static int MIDDLE = 2;
        public static int END = 3;
    }
}
