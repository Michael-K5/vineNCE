#' Import packages
#' @name imports
#' @keywords internal
#' @import rvinecopulib
#' @import tensorflow
#' @import keras
#' @importFrom graphics box contour hist mtext par
#' @importFrom stats cor poly predict qnorm runif sd
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

# register the binary cross entropy metric on library loading
.onLoad <- function(libname, pkgname) {
  reticulate::py_run_string("
import tensorflow as tf
from tensorflow import keras

@keras.saving.register_keras_serializable(name='binary_accuracy_from_logits')
def binary_accuracy_from_logits(y_true, y_pred):
    y_pred = tf.nn.sigmoid(y_pred)
    return keras.metrics.binary_accuracy(y_true, y_pred)
")
  assign("binary_accuracy_from_logits", reticulate::py$binary_accuracy_from_logits,
         envir = parent.env(environment()))
}
