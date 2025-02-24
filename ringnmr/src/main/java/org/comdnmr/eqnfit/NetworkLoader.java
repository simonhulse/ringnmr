package org.comdnmr.eqnfit;

import java.io.File;
import java.io.IOException;

import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

import java.util.HashMap;
import java.util.Map;

import org.tensorflow.SavedModelBundle;

import org.nmrfx.utilities.UnZipper;

/**
 *
 * @author simonhulse
 */

// TODO: currently, NetworkLoader is specifically used to load CPMG parameter
// prediction NNs, however it could be generalised to any zipfile containing
// tensorflow SavedModelBundle models.
class NetworkLoader {
    private static NetworkLoader instance = null;

    private Map<String, SavedModelBundle> loadedNetworks;

    private String name = "parameter-networks";

    private File zipSrc = new File(String.format("/data/%s.zip", name));
    private File zipDst = new File(System.getProperty("java.io.tmpdir"), String.format("%s.zip", name));
    private File zipDstParent = zipDst.getParentFile();
    private File rootDir = new File(zipDstParent, name);

    private NetworkLoader() throws IOException {
        if (!rootDir.exists()) {
            // TODO: Gracefully handle IOException
            Files.copy(
                getClass().getResourceAsStream(zipSrc.getPath()),
                zipDst.toPath(),
                StandardCopyOption.REPLACE_EXISTING);

            UnZipper unZipper = new UnZipper(zipDstParent, zipDst.getPath());
            unZipper.unzip();

            Files.deleteIfExists(zipDst.toPath());
        }

        loadedNetworks = new HashMap<>();
    }

    public static synchronized NetworkLoader getNetworkLoader() throws IOException {
        if (instance == null) {
            instance = new NetworkLoader();
        }

        return instance;
    }

    public SavedModelBundle fetchNetwork(String path) {
        SavedModelBundle network;
        if (loadedNetworks.containsKey(path)) {
            network = loadedNetworks.get(path);
        } else {
            File absolutePath = getAbsolutePath(path);
            network = SavedModelBundle.load(absolutePath.getPath(), "serve");
            loadedNetworks.put(path, network);
        }
        return network;
    }

    private File getAbsolutePath(String relativePath) {
        return new File(rootDir, relativePath);
    }
}
