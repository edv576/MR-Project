using System;
using System.IO;
using System.Linq;
using Assimp;

namespace MRAnalize
{
    class Program
    {
        static void Main(string[] args)
        {
            var path = Directory.GetCurrentDirectory();
            var files = Directory
                //get all directories
                .GetDirectories(Path.Combine(path, "PSB"))
                //and all files within
                .Select(dir => (category: new DirectoryInfo(dir).Name, files: Directory.GetFiles(dir)))
                //flatten the files and combine with the category
                .SelectMany(dirInfo => dirInfo.files.Select(file => (dirInfo.category, file)))
                //select only the .off files
                .Where(fileInfo => Path.GetExtension(fileInfo.file) == ".off");

            //open file to write to
            using var database = new StreamWriter(Path.Combine(path, $"out/db{DateTime.Now.Hour}{DateTime.Now.Minute}{DateTime.Now.Second}.csv"));

            //write header
            database.WriteLine("filename,class,vertices,faces");

            foreach (var (category, file) in files)
            {
                //import assets
                var importContext = new AssimpContext();
                //triangluate them
                var scene = importContext.ImportFile(file);
                //foreach mesh get the total number of vertices and faces and sum them
                var nrvertices = scene.Meshes.Select(m => m.Vertices.Count).Sum();
                var nrfaces = scene.Meshes.Select(m => m.Faces.Count).Sum();

                //add the to the database
                database.WriteLine($"{Path.GetFileName(file)},{category},{nrvertices},{nrfaces}");

                //safe the triangulated file
                if (!importContext.ExportFile(scene, Path.Combine(path, $"out/{Path.GetFileNameWithoutExtension(file)}.ply"), "ply"))
                    Console.WriteLine("Export failed");
            }
        }
    }
}
