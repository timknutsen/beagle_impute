/* Sandve Lab — main.js */

(function () {
  'use strict';

  const toggle   = document.getElementById('navToggle');
  const menu     = document.getElementById('mobileMenu');

  if (!toggle || !menu) return;

  toggle.addEventListener('click', function () {
    const isOpen = menu.classList.toggle('open');
    toggle.setAttribute('aria-expanded', isOpen);

    // Animate burger → X
    const spans = toggle.querySelectorAll('span');
    if (isOpen) {
      spans[0].style.transform = 'translateY(7px) rotate(45deg)';
      spans[1].style.opacity   = '0';
      spans[2].style.transform = 'translateY(-7px) rotate(-45deg)';
    } else {
      spans[0].style.transform = '';
      spans[1].style.opacity   = '';
      spans[2].style.transform = '';
    }
  });

  // Close menu when a link is clicked
  menu.querySelectorAll('.nav__mobile-link').forEach(function (link) {
    link.addEventListener('click', function () {
      menu.classList.remove('open');
      toggle.setAttribute('aria-expanded', 'false');
      const spans = toggle.querySelectorAll('span');
      spans[0].style.transform = '';
      spans[1].style.opacity   = '';
      spans[2].style.transform = '';
    });
  });

  // Scroll-based nav shadow
  window.addEventListener('scroll', function () {
    document.querySelector('.nav').style.boxShadow =
      window.scrollY > 10 ? '0 2px 20px rgba(0,0,0,.25)' : '';
  }, { passive: true });
}());
