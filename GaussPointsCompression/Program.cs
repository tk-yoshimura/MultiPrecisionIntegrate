using System.IO.Compression;

const string dirpath_src = "../../../../gauss_points/";
const string dirpath_dst = "../../../../MultiPrecisionIntegrate/Resource/";

foreach (string filepath in Directory.EnumerateFiles(dirpath_src, "*.bin")) {
    string filename = Path.GetFileName(filepath);

    using FileStream sfs = File.OpenRead(filepath);
    using FileStream dfs = File.Create(dirpath_dst + filename);

    byte[] data = new byte[sfs.Length];
    sfs.Read(data, 0, checked((int)sfs.Length));

    using var compressor = new GZipStream(dfs, CompressionLevel.SmallestSize);

    compressor.Write(data, 0, data.Length);
}