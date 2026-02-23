def conf_mfx(coords_input, conf_offset, params, *args, **kwargs):
        "Pre-calibrated coordinate transform for etMINFLUX."

        conf_size_x_um = params[0]
        conf_size_y_um = params[1]
        conf_size_x_px = params[2]
        conf_size_y_px = params[3]

        x_um = coords_input[0]*conf_size_x_um/conf_size_x_px
        y_um = coords_input[1]*conf_size_y_um/conf_size_y_px
        x_um_corr = x_um - conf_size_x_um/2  # local in image
        y_um_corr = y_um - conf_size_y_um/2  # local in image
        x_um_glob = x_um_corr + conf_offset[0]  # global position, um
        y_um_glob = y_um_corr + conf_offset[1]  # global position, um

        return [x_um_glob, y_um_glob]