namespace MultiPrecisionIntegrate {
    [global::System.CodeDom.Compiler.GeneratedCodeAttribute("System.Resources.Tools.StronglyTypedResourceBuilder", "17.0.0.0")]
    [global::System.Diagnostics.DebuggerNonUserCodeAttribute()]
    [global::System.Runtime.CompilerServices.CompilerGeneratedAttribute()]
    internal class Resource {
        
        private static global::System.Resources.ResourceManager resourceMan;
        
        private static global::System.Globalization.CultureInfo resourceCulture;
        
        [global::System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1811:AvoidUncalledPrivateCode")]
        internal Resource() {
        }

        [global::System.ComponentModel.EditorBrowsableAttribute(global::System.ComponentModel.EditorBrowsableState.Advanced)]
        internal static global::System.Resources.ResourceManager ResourceManager {
            get {
                if (object.ReferenceEquals(resourceMan, null)) {
                    global::System.Resources.ResourceManager temp = new global::System.Resources.ResourceManager("MultiPrecisionIntegrate.Resource", typeof(Resource).Assembly);
                    resourceMan = temp;
                }
                return resourceMan;
            }
        }
        
        [global::System.ComponentModel.EditorBrowsableAttribute(global::System.ComponentModel.EditorBrowsableState.Advanced)]
        internal static global::System.Globalization.CultureInfo Culture {
            get {
                return resourceCulture;
            }
            set {
                resourceCulture = value;
            }
        }
        
        internal static byte[] hermite_n64 {
            get {
                object obj = ResourceManager.GetObject("hermite_n64", resourceCulture);
                return ((byte[])(obj));
            }
        }
        
        internal static byte[] laguerre_n64 {
            get {
                object obj = ResourceManager.GetObject("laguerre_n64", resourceCulture);
                return ((byte[])(obj));
            }
        }
        
        internal static byte[] legendre_01_n64 {
            get {
                object obj = ResourceManager.GetObject("legendre_01_n64", resourceCulture);
                return ((byte[])(obj));
            }
        }
    }
}
