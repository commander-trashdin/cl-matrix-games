;;;; matrix-reloaded.lisp

(in-package #:matrix-reloaded)

(defun equidimensional (a)
  (or (< (array-rank a) 2)
      (apply #'= (array-dimensions a))))

(deftype matrix (&optional type horsize versize)
  `(array ,type ,(list horsize versize)))

(deftype square-matrix (&optional type size)
  `(matrix ,type ,size ,size))

(defun lambda-e-matrix (size &optional (coef 1))
  (let ((lem (make-array (list size size) :initial-element 0 :element-type 'fixnum)))
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

(declaim (ftype (function (matrix matrix) matrix) matprod))
(defun matprod (fst snd)
  (destructuring-bind (m n) (array-dimensions fst)
    (destructuring-bind (n1 l) (array-dimensions snd)
      (assert (and (= n n1)
                   (eql (array-element-type fst) (array-element-type fst))))
      (let ((prod (make-array (list m l) :element-type (array-element-type fst))))
          (loop :for i :from 0 :below m :do
            (loop :for j :from 0 :below l :do
              (setf (aref prod i j)
                    (loop :for k :from 0 :below n
                          :sum (* (aref fst i k) (aref snd k j))))))
        prod))))


(defun matexpt (mat deg)
  (cond
    ((zerop deg) (lambda-e-matrix (array-dimension mat 1)))
    ((= 1 deg) mat)
    (t (let ((prevdeg (matexpt mat (ash deg -1))))
         (if (evenp deg)
             (matprod prevdeg prevdeg)
             (matprod mat (matprod prevdeg prevdeg)))))))

(defun %base-vect (dim pos)
  (let ((res (make-array dim :element-type 'fixnum :initial-element 0)))
    (setf (aref res pos) 1)
    res))

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
