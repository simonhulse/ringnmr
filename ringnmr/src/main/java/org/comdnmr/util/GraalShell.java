package org.comdnmr.util;

import org.graalvm.polyglot.*;

import java.io.*;
import java.util.Set;

public class GraalShell {

    public static void main(String[] argv) {
        System.out.println("argv " + argv.length);

        if (argv.length > 0) {
            Context context = Context.newBuilder().allowAllAccess(true)
                    .option("python.EmulateJython", "true").build();
            var pyLang = context.getEngine().getLanguages().get("python");

            for (var option : pyLang.getOptions()) {
                System.out.println("opt " + option.getName() + " " + option.getUsageSyntax());
            }
            File file = new File(argv[0]);
            Source source = Source.newBuilder("python", file)
                    .buildLiteral();
            context.eval(source);
        } else {
            BufferedReader input = new BufferedReader(new InputStreamReader(System.in));
            PrintStream output = System.out;
            Context context = Context.newBuilder().allowAllAccess(true).build();
            Set<String> languages = context.getEngine().getLanguages().keySet();
            output.println("Shell for " + languages + ":");
            String language = languages.iterator().next();
            for (; ; ) {
                try {
                    output.print(language + "> ");
                    String line = input.readLine();
                    if (line == null) {
                        break;
                    } else if (languages.contains(line)) {
                        language = line;
                        continue;
                    }
                    Source source = Source.newBuilder("python", line, "a")
                            .buildLiteral();
                    context.eval(source);
                } catch (PolyglotException t) {
                    if (t.isExit()) {
                        break;
                    }
                    t.printStackTrace();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }

    }

}
