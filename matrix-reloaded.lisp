;;;; matrix-reloaded.lisp

(in-package #:matrix-reloaded)

#+sbcl (require :sb-gmp)

(declaim (optimize (speed 3) (debug 0) (space 0) (safety 0)))

(defun equidimensional (a)
  (or (< (array-rank a) 2)
      (apply #'= (array-dimensions a))))

(deftype matrix (&optional type horsize versize)
  `(array ,type ,(list horsize versize)))

(deftype square-matrix (&optional type size)
  `(matrix ,type ,size ,size))

(defun lambda-e-matrix (size &optional (coef 1))
  (let ((lem (make-array (list size size) :initial-element 0 :element-type 'number)))
    (loop :for i :from 0 :below size
          :do (setf (aref lem i i) coef))
    lem))

(declaim (ftype (function (matrix matrix) matrix) matsum))
(defun matsum (fst snd)
  (assert (and (equalp (array-dimensions fst) (array-dimensions snd))
               (eql (array-element-type fst) (array-element-type fst))))
  (let ((sum (make-array (array-dimensions fst) :element-type (array-element-type fst))))
    (destructuring-bind (n m) (array-dimensions sum)
      (loop :for i :from 0 :below n :do
        (loop :for j :from 0 :below m do
          (setf (aref sum i j) (+ (aref fst i j) (aref snd i j))))))
    sum))

(declaim (ftype (function (matrix matrix &optional (or null fixnum) (or null fixnum) (or null fixnum)) matrix) matprod))
(defun matprod (fst snd &optional d1 d2 d3)
  (declare (optimize (speed 3) (debug 0) (space 0) (safety 0)))
  (if d1
     (let ((prod (make-array (list d1 d3) :element-type (array-element-type fst))))
        (loop :for i :from 0 :below d1 :do
          (loop :for j :from 0 :below d3 :do
            (setf (aref prod i j)
                  (loop :for k :from 0 :below d2
                        :sum (the integer (* (the integer (aref fst i k)) (the integer (aref snd k j))))))))
      prod)
     (destructuring-bind (m n) (array-dimensions fst)
       (destructuring-bind (n1 l) (array-dimensions snd)
         (assert (and (= n n1)
                      (eql (array-element-type fst) (array-element-type fst))))
         (let ((prod (make-array (list m l) :element-type (array-element-type fst))))
             (loop :for i :from 0 :below m :do
               (loop :for j :from 0 :below l :do
                 (setf (aref prod i j)
                       (loop :for k :from 0 :below n
                             :sum (the integer (* (the integer (aref fst i k)) (the integer (aref snd k j))))))))
           prod)))))

(declaim (ftype (function (matrix fixnum &optional fixnum) matrix) matexpt))
(defun matexpt (mat deg &optional dim)
  (declare (optimize (speed 3) (debug 0) (space 0) (safety 0)))
  (cond
    ((zerop deg) (lambda-e-matrix (array-dimension mat 1)))
    ((= 1 deg) mat)
    (t (let* ((d (or dim (array-dimension mat 1)))
              (prevdeg (the matrix (matexpt mat (ash deg -1) d))))
         (if (logbitp 0 deg)
             (matprod mat (matprod prevdeg prevdeg d d d) d d d)
             (matprod prevdeg prevdeg d d d))))))

(declaim (ftype (function (matrix fixnum) matrix) matexpt1))
(defun matexpt1 (mat deg)
  (do ((res (lambda-e-matrix (array-dimension mat 1))))
      ((< deg 1) res)
      (when (logbitp 0 deg)
        (setf res (matprod res mat)))
      (setf mat (the matrix (matprod mat mat))
            deg (the fixnum (ash deg -1)))))


(declaim (ftype (function (fixnum fixnum) (simple-vector *)) %base-vect)
         (inline %base-vect))
(defun %base-vect (dim pos)
  (let ((res (make-array dim :element-type 'fixnum :initial-element 0)))
    (setf (aref res pos) 1)
    res))


(declaim (ftype (function (fixnum fixnum) list) %base-list)
         (inline %base-list))
(defun %base-list (dim pos)
  (let ((res (make-list dim :initial-element 0)))
    (setf (elt res pos) 1)
    res))

(defun solve-reccur (index rule &rest base)
  (let* ((l (length base))
         (matbase (make-array (list l l) :element-type 'integer :initial-element 0))
         (matdenom (make-array (list l l) :element-type 'integer :initial-element 0)))
    (loop :for i :from 0 :below l
           :do (setf (aref matbase 0 i) (elt base i)
                     (aref matdenom i (1- l)) (apply rule (%base-list l i))))
    (loop :for j :from 0 :below (1- l)
          :do (setf (aref matdenom (1+ j) j) 1))
    (if (<= index l)
        matbase
        (matprod matbase (matexpt matdenom (- index l))))))
