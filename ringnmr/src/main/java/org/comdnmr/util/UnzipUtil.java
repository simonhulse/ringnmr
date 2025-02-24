package org.comdnmr.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class UnzipUtil {

    public static void unzip(String srcPath, String dstDir) {
        try (ZipInputStream zis = new ZipInputStream(new FileInputStream(srcPath))) {
            ZipEntry zipEntry = zis.getNextEntry();

            while (zipEntry != null) {
                String fileName = zipEntry.getName();
                String filePath = dstDir + "/" + fileName;

                if (!zipEntry.isDirectory()) {
                    // Create parent directories if necessary
                    new File(filePath).getParentFile().mkdirs();

                    try (FileOutputStream fos = new FileOutputStream(filePath)) {
                        byte[] buffer = new byte[1024];
                        int len;
                        while ((len = zis.read(buffer)) > 0) {
                            fos.write(buffer, 0, len);
                        }
                    }
                }

                zipEntry = zis.getNextEntry();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
