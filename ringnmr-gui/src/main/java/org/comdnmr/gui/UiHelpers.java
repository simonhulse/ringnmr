package org.comdnmr.gui;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;

import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Label;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.util.Duration;

public class UiHelpers {

    private static final int DEFAULT_BOX_GAP = 2;
    private static final int DEFAULT_PANE_GAP = 10;
    private static final Pos DEFAULT_ALIGNMENT = Pos.CENTER_LEFT;
    private static final int DEFAULT_SPACER_WIDTH = 5;

    private static final String HELPER_MESSAGES_DIR = "helper_messages";

    public static GridPane createDefaultGridPane() {
        GridPane pane = new GridPane();
        pane.setHgap(DEFAULT_PANE_GAP);
        pane.setVgap(DEFAULT_PANE_GAP);
        pane.setAlignment(DEFAULT_ALIGNMENT);
        return pane;
    }

    public static Region createSpacer() {
        return createSpacer(DEFAULT_SPACER_WIDTH);
    }

    public static Region createSpacer(int width) {
        Region spacer = new Region();
        spacer.setPrefWidth(width);
        return spacer;
    }

    public static HBox createDefaultHBox(Node... nodes) {
        HBox hBox = new HBox(DEFAULT_BOX_GAP, nodes);
        hBox.setAlignment(DEFAULT_ALIGNMENT);
        return hBox;
    }

    /** Create a small helper icon that shows a message on hover/focus. */
    public static Label createHelperIcon(String fileName) {
        Label icon = new Label("?");
        icon.getStyleClass().add("helper-icon");

        String message = getMessage(fileName);
        Tooltip tooltip = createTooltip(message);
        Tooltip.install(icon, tooltip);

        // keyboard accessible: show tooltip when focused
        icon.focusedProperty().addListener((obs, oldV, newV) -> {
            if (newV) tooltip.show(
                icon,
                icon
                    .localToScreen(icon.getBoundsInLocal())
                    .getMinX(),
                icon
                    .localToScreen(icon.getBoundsInLocal())
                    .getMaxY()
            );
            else tooltip.hide();
        });

        return icon;
    }

    public static HBox createElementWithHelper(Node element, String fileName) {
        HBox hBox = new HBox();
        hBox.setSpacing(5);
        hBox.setAlignment(Pos.CENTER_LEFT);
        Label helper = createHelperIcon(fileName);
        hBox.getChildren().addAll(element, helper);
        return hBox;
    }

    private static String getMessage(String fileName) {
        String path = "/" + HELPER_MESSAGES_DIR + "/" + fileName;
        try (InputStream in = UiHelpers.class.getResourceAsStream(path)) {
            if (in == null) throw new IllegalArgumentException("Resource not found: " + path);
            return new String(in.readAllBytes(), StandardCharsets.UTF_8);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static Tooltip createTooltip(String message) {
        Tooltip tooltip = new Tooltip(message);
        tooltip.setShowDelay(Duration.ZERO);
        tooltip.setShowDuration(Duration.INDEFINITE);
        tooltip.setHideDelay(Duration.ZERO);
        tooltip.setWrapText(true);
        tooltip.setMaxWidth(300);
        return tooltip;
    }
}
